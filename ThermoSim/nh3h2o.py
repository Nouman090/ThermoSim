"""
nh3h2o.py -- CoolProp-style thermodynamic property library for NH3-H2O mixtures.

Give any valid pair of properties plus the overall NH3 mass fraction; get the
complete state back.

    >>> import nh3h2o as aw
    >>> aw.PropsSI('T', 'P', 10, 'Q', 0, 0.4)      # bubble temperature
    >>> st = aw.state(z=0.4, P=10, H=500)          # full state dict
    >>> st.T, st.Q, st.x, st.y

MODEL
  Single phase   Ziegler & Trepp (1984) Gibbs free-energy formulation.
  Equilibrium    solved from the SAME Gibbs function (equal chemical
                 potentials), so enthalpy, entropy and VLE are mutually
                 consistent and energy balances close exactly.
  Optional       Patek & Klomfar (1995) explicit bubble/dew correlations,
                 selectable with vle='PK'. Closer to experimental VLE data but
                 NOT thermodynamically consistent with the enthalpies.

UNITS   T [K], P [bar], h [kJ/kg], s [kJ/kg/K], v [m3/kg], u [kJ/kg],
        cp [kJ/kg/K], all compositions are NH3 MASS fractions.
REFERENCE   h = s = 0 for the saturated liquid of each pure component at 273.15 K.

VALIDITY    230-500 K, 0.2-50 bar (Ziegler & Trepp validated envelope).
            Accuracy degrades rapidly above 50 bar; see check_range().

NOTE    Vendored into ThermoSim (units here are still bar/kJ, as above; the
        `Prop`/`NH3H2OFluid` layer in state.py converts to/from ThermoSim's
        SI convention). This module's own units and public API are otherwise
        unchanged from the standalone library.
"""
from __future__ import annotations
import warnings
from dataclasses import dataclass, asdict

import numpy as np
import sympy as sp
from scipy.optimize import brentq

__all__ = ['state', 'PropsSI', 'bubble_T', 'dew_T', 'bubble_P', 'dew_P',
           'props_single_phase', 'AmmoniaWaterError', 'T_MIN', 'T_MAX', 'P_MIN', 'P_MAX',
           'NH3H2OFluid']


class AmmoniaWaterError(RuntimeError):
    """Raised when a state cannot be resolved."""


# --------------------------------------------------------------------------
# constants and Ziegler-Trepp coefficients
# --------------------------------------------------------------------------
M_NH3, M_H2O = 17.03, 18.015          # kg/kmol
R, TB, PB = 8.314, 100.0, 10.0        # kJ/kmol/K, K, bar
T_MIN, T_MAX, P_MIN, P_MAX = 230.0, 500.0, 0.2, 50.0
_TINY = 1e-12

_CO = {
    'ammonia': dict(A1=3.971423e-2, A2=-1.790557e-5, A3=-1.308905e-2, A4=3.752836e-3,
                    B1=1.634519e1, B2=-6.508119, B3=1.448937,
                    C1=-1.049377e-2, C2=-8.288224, C3=-6.647257e2, C4=-3.045352e3,
                    D1=3.673647, D2=9.989629e-2, D3=3.617622e-2,
                    hL0=4.878573, hg0=26.468873, sL0=1.644773, sg0=8.339026,
                    T0=3.2252, P0=2.0),
    'water':   dict(A1=2.748796e-2, A2=-1.016665e-5, A3=-4.452025e-3, A4=8.389246e-4,
                    B1=1.214557e1, B2=-1.898065, B3=2.911966e-1,
                    C1=2.136131e-2, C2=-3.169291e1, C3=-4.634611e4, C4=0.0,
                    D1=4.019170, D2=-5.175550e-2, D3=1.951939e-2,
                    hL0=21.821141, hg0=60.965058, sL0=5.733498, sg0=13.453430,
                    T0=5.0705, P0=3.0),
}
_E = [-41.733398, 0.02414, 6.702285, -0.011475, 63.608967, -62.490768,
      1.761064, 0.008626, 0.387983, -0.004772, -4.648107, 0.836376,
      -3.553627, 0.000904, 24.361723, -20.736547]


# --------------------------------------------------------------------------
# symbolic setup -- differentiated and lambdified ONCE at import
# --------------------------------------------------------------------------
_Tr, _Pr, _xs = sp.symbols('Tr Pr xs', positive=True)


def _pure_gibbs(fluid):
    c = _CO[fluid]
    T0, P0 = c['T0'], c['P0']
    gL = (c['hL0'] - _Tr * c['sL0'] + c['B1'] * (_Tr - T0)
          + c['B2'] / 2 * (_Tr**2 - T0**2) + c['B3'] / 3 * (_Tr**3 - T0**3)
          - c['B1'] * _Tr * sp.log(_Tr / T0) - c['B2'] * _Tr * (_Tr - T0)
          - c['B3'] / 2 * _Tr * (_Tr**2 - T0**2)
          + (c['A1'] + c['A3'] * _Tr + c['A4'] * _Tr**2) * (_Pr - P0)
          + c['A2'] / 2 * (_Pr**2 - P0**2))
    gG = (c['hg0'] - _Tr * c['sg0'] + c['D1'] * (_Tr - T0)
          + c['D2'] / 2 * (_Tr**2 - T0**2) + c['D3'] / 3 * (_Tr**3 - T0**3)
          - c['D1'] * _Tr * sp.log(_Tr / T0) - c['D2'] * _Tr * (_Tr - T0)
          - c['D3'] / 2 * _Tr * (_Tr**2 - T0**2)
          + _Tr * sp.log(_Pr / P0) + c['C1'] * (_Pr - P0)
          + c['C2'] * (_Pr / _Tr**3 - 4 * P0 / T0**3 + 3 * P0 * _Tr / T0**4)
          + c['C3'] * (_Pr / _Tr**11 - 12 * P0 / T0**11 + 11 * P0 * _Tr / T0**12)
          + c['C4'] / 3 * (_Pr**3 / _Tr**11 - 12 * P0**3 / T0**11
                           + 11 * P0**3 * _Tr / T0**12))
    return gL, gG


_GLa, _GGa = _pure_gibbs('ammonia')
_GLw, _GGw = _pure_gibbs('water')

_F1 = _E[0] + _E[1] * _Pr + (_E[2] + _E[3] * _Pr) * _Tr + _E[4] / _Tr + _E[5] / _Tr**2
_F2 = _E[6] + _E[7] * _Pr + (_E[8] + _E[9] * _Pr) * _Tr + _E[10] / _Tr + _E[11] / _Tr**2
_F3 = _E[12] + _E[13] * _Pr + _E[14] / _Tr + _E[15] / _Tr**2
_GE = (_F1 + _F2 * (2 * _xs - 1) + _F3 * (2 * _xs - 1)**2) * (1 - _xs) * _xs

_L = lambda e: sp.lambdify((_Tr, _Pr, _xs), e, 'numpy')
_f = {}
for _n, _G in [('La', _GLa), ('Ga', _GGa), ('Lw', _GLw), ('Gw', _GGw)]:
    _f[_n + '_g'] = _L(_G)
    _f[_n + '_h'] = _L(sp.diff(_G / _Tr, _Tr))
    _f[_n + '_s'] = _L(sp.diff(_G, _Tr))
    _f[_n + '_v'] = _L(sp.diff(_G, _Pr))
_f['E_g'] = _L(_GE)
_f['E_h'] = _L(sp.diff(_GE / _Tr, _Tr))
_f['E_s'] = _L(sp.diff(_GE, _Tr))
_f['E_v'] = _L(sp.diff(_GE, _Pr))
_f['E_dx'] = _L(sp.diff(_GE, _xs))          # analytic dG_E/dx -> exact chem. potentials


# --------------------------------------------------------------------------
# composition helpers
# --------------------------------------------------------------------------
def to_mole(X):
    """NH3 mass fraction -> NH3 mole fraction."""
    X = min(max(float(X), 0.0), 1.0)
    return (X / M_NH3) / (X / M_NH3 + (1.0 - X) / M_H2O)


def to_mass(x):
    """NH3 mole fraction -> NH3 mass fraction."""
    x = min(max(float(x), 0.0), 1.0)
    return x * M_NH3 / (x * M_NH3 + (1.0 - x) * M_H2O)


def _Mbar(x):
    return x * M_NH3 + (1.0 - x) * M_H2O


def _xlogx(x):
    return 0.0 if (x <= 0.0 or x >= 1.0) else x * np.log(x)


def _s_ideal(x):
    """Ideal mixing entropy, kJ/kmol/K. NaN-safe at the pure limits."""
    return -R * (_xlogx(x) + _xlogx(1.0 - x))


def check_range(T, P, warn=True):
    ok = (T_MIN <= T <= T_MAX) and (P_MIN <= P <= P_MAX)
    if warn and not ok:
        warnings.warn(
            f'T={T:.1f} K, P={P:.2f} bar lies outside the Ziegler-Trepp validated '
            f'envelope ({T_MIN}-{T_MAX} K, {P_MIN}-{P_MAX} bar). Vapour-phase '
            f'enthalpy error grows rapidly above 50 bar.',
            RuntimeWarning, stacklevel=3)
    return ok


# --------------------------------------------------------------------------
# single-phase properties
# --------------------------------------------------------------------------
def props_single_phase(T, P, X, phase):
    """Specific properties of one homogeneous phase.

    phase : 'liquid' or 'vapor'; X is the NH3 mass fraction OF THAT PHASE.
    Returns dict with h, s, v, u [kJ/kg, kJ/kg/K, m3/kg, kJ/kg].
    """
    x = to_mole(X)
    tr, pr = T / TB, P / PB
    M = _Mbar(x)
    sid = _s_ideal(x)

    if phase.startswith('l'):
        ha = -R * TB * tr**2 * _f['La_h'](tr, pr, x)
        hw = -R * TB * tr**2 * _f['Lw_h'](tr, pr, x)
        sa = -R * _f['La_s'](tr, pr, x)
        sw = -R * _f['Lw_s'](tr, pr, x)
        va = R * TB / PB * _f['La_v'](tr, pr, x)
        vw = R * TB / PB * _f['Lw_v'](tr, pr, x)
        hE = -R * TB * tr**2 * _f['E_h'](tr, pr, x)
        sE = -R * _f['E_s'](tr, pr, x)
        vE = R * TB / PB * _f['E_v'](tr, pr, x)
    else:
        ha = -R * TB * tr**2 * _f['Ga_h'](tr, pr, x)
        hw = -R * TB * tr**2 * _f['Gw_h'](tr, pr, x)
        sa = -R * _f['Ga_s'](tr, pr, x)
        sw = -R * _f['Gw_s'](tr, pr, x)
        va = R * TB / PB * _f['Ga_v'](tr, pr, x)
        vw = R * TB / PB * _f['Gw_v'](tr, pr, x)
        hE = sE = vE = 0.0                      # vapour treated as an ideal mixture

    h = (x * ha + (1 - x) * hw + hE) / M
    s = (x * sa + (1 - x) * sw + sid + sE) / M
    v = (x * va + (1 - x) * vw + vE) / M * 0.01  # kJ/(kmol.bar) -> m3/kmol
    return dict(h=h, s=s, v=v, u=h - P * v * 100.0)


def _g_pure(T, P, which):
    """Pure-component molar Gibbs energy [kJ/kmol]. which in {La,Ga,Lw,Gw}."""
    return R * TB * _f[which + '_g'](T / TB, P / PB, 0.5)


# --------------------------------------------------------------------------
# vapour-liquid equilibrium
# --------------------------------------------------------------------------
def _lnK(T, P, x):
    """Natural logs of the K-values at liquid mole fraction x.

    Equal chemical potentials give
        mu_A = g_A^L + RT ln x + G_E + (1-x) dG_E/dx  =  g_A^V + RT ln y
        mu_W = g_W^L + RT ln(1-x) + G_E - x dG_E/dx   =  g_W^V + RT ln(1-y)
    so K_A = y/x and K_W = (1-y)/(1-x) depend on (T, P, x) ONLY -- not on y.
    That is what turns every flash into a sequence of 1-D root finds.
    """
    tr, pr = T / TB, P / PB
    RT = R * T
    GE = R * TB * _f['E_g'](tr, pr, x)
    dGE = R * TB * _f['E_dx'](tr, pr, x)
    lnKA = (_g_pure(T, P, 'La') - _g_pure(T, P, 'Ga') + GE + (1 - x) * dGE) / RT
    lnKW = (_g_pure(T, P, 'Lw') - _g_pure(T, P, 'Gw') + GE - x * dGE) / RT
    return lnKA, lnKW


def _bubble_residual(T, P, x):
    """Zero when (T, P, x) sits on the bubble surface: sum of vapour fractions = 1."""
    lnKA, lnKW = _lnK(T, P, x)
    return np.exp(lnKA) * x + np.exp(lnKW) * (1 - x) - 1.0


def _y_from_x(T, P, x):
    lnKA, _ = _lnK(T, P, x)
    return min(max(np.exp(lnKA) * x, 0.0), 1.0)


def _solve_bracketed(fun, lo, hi, n=48, name='root', log=False, pick='first'):
    """Scan for a sign change then bisect. Slower than a Newton step but it
    does not diverge, which matters more inside nested flash loops.

    pick='last' returns the highest root. The Ziegler-Trepp vapour Gibbs
    function contains C3/Tr**11 and C4*Pr**3/Tr**11 terms that blow up at low
    reduced temperature and can manufacture spurious roots there, so the
    saturation solvers take the highest (physical) crossing.
    """
    grid = (np.geomspace(max(lo, 1e-6), hi, n) if log else np.linspace(lo, hi, n))
    prev_t, prev_v = None, None
    brackets = []
    for t in grid:
        try:
            v = fun(t)
        except (ValueError, FloatingPointError, OverflowError, AmmoniaWaterError):
            prev_t, prev_v = None, None
            continue
        if not np.isfinite(v):
            prev_t, prev_v = None, None
            continue
        if prev_v is not None and np.sign(v) != np.sign(prev_v):
            if pick == 'first':
                return brentq(fun, prev_t, t, xtol=1e-10, rtol=1e-12, maxiter=200)
            brackets.append((prev_t, t))
        prev_t, prev_v = t, v
    if brackets:
        a, b = brackets[-1]
        return brentq(fun, a, b, xtol=1e-10, rtol=1e-12, maxiter=200)
    raise AmmoniaWaterError(f'no {name} bracketed in [{lo:.4g}, {hi:.4g}]')


def _all_roots(fun, lo, hi, n=110, log=False):
    """Every sign-change root in the interval, plus the observed span of the
    residual. A tiny span means the two inputs do not actually determine the
    state (e.g. compressed-liquid enthalpy is nearly independent of pressure)."""
    grid = (np.geomspace(max(lo, 1e-6), hi, n) if log else np.linspace(lo, hi, n))
    roots, vals = [], []
    prev_t, prev_v = None, None
    for t in grid:
        try:
            v = fun(t)
        except (ValueError, FloatingPointError, OverflowError, AmmoniaWaterError):
            prev_t, prev_v = None, None
            continue
        if not np.isfinite(v):
            prev_t, prev_v = None, None
            continue
        vals.append(v)
        if prev_v is not None and np.sign(v) != np.sign(prev_v):
            roots.append(brentq(fun, prev_t, t, xtol=1e-10, rtol=1e-12, maxiter=200))
        prev_t, prev_v = t, v
    span = (max(vals) - min(vals)) if vals else 0.0
    return roots, span


def bubble_T(P, z, vle='ZT'):
    """Bubble-point temperature [K] at pressure P and liquid NH3 mass fraction z."""
    if vle == 'PK':
        return _pk_bubble_T(P, z)
    x = to_mole(z)
    if x <= _TINY or x >= 1 - _TINY:                      # pure component
        return _pure_Tsat(P, 'ammonia' if x > 0.5 else 'water')
    Ta, Tw = _envelope(P)
    return _solve_bracketed(lambda T: _bubble_residual(T, P, x),
                            Ta + 1e-6, Tw - 1e-6, name='bubble temperature')


def _x_at_TP(T, P):
    """Liquid NH3 mole fraction on the bubble surface at (T, P)."""
    Ta, Tw = _envelope(P)
    if not (Ta < T < Tw):
        raise AmmoniaWaterError(f'no two-phase region at T={T:.3f} K, P={P:.4g} bar')
    f = lambda x: _bubble_residual(T, P, x)
    lo, hi = 1e-9, 1 - 1e-9
    try:
        # between the pure saturation temperatures the residual is negative at the
        # water end and positive at the ammonia end, so a bracket always exists
        if f(lo) * f(hi) < 0:
            return brentq(f, lo, hi, xtol=1e-12, rtol=1e-13, maxiter=100)
    except (ValueError, FloatingPointError, OverflowError):
        pass
    return _solve_bracketed(f, lo, hi, n=64, name='bubble composition')


def dew_T(P, z, vle='ZT'):
    """Dew-point temperature [K] at pressure P and vapour NH3 mass fraction z."""
    if vle == 'PK':
        return _pk_dew_T(P, z)
    y_t = to_mole(z)
    if y_t <= _TINY or y_t >= 1 - _TINY:
        return _pure_Tsat(P, 'ammonia' if y_t > 0.5 else 'water')

    Ta, Tw = _envelope(P)

    def resid(T):
        return _y_from_x(T, P, _x_at_TP(T, P)) - y_t

    return _solve_bracketed(resid, Ta + 1e-6, Tw - 1e-6, n=48,
                            name='dew temperature')


def bubble_P(T, z, vle='ZT'):
    return _solve_bracketed(lambda P: bubble_T(P, z, vle) - T, 0.02, 200.0,
                            n=60, name='bubble pressure', log=True)


def dew_P(T, z, vle='ZT'):
    return _solve_bracketed(lambda P: dew_T(P, z, vle) - T, 0.02, 200.0,
                            n=60, name='dew pressure', log=True)


_TSAT_CACHE = {}


def _pure_Tsat(P, fluid):
    """Saturation temperature of a pure component from g_L = g_V. Cached."""
    key = (round(float(P), 10), fluid)
    if key in _TSAT_CACHE:
        return _TSAT_CACHE[key]
    a = 'a' if fluid == 'ammonia' else 'w'
    T = _solve_bracketed(
        lambda T: _g_pure(T, P, 'L' + a) - _g_pure(T, P, 'G' + a),
        150.0, 700.0, n=140, name=f'{fluid} saturation temperature', pick='last')
    _TSAT_CACHE[key] = T
    return T


def _envelope(P):
    """(T_sat NH3, T_sat H2O) at pressure P -- the two-phase region lies between
    them for every intermediate composition, which gives every temperature root
    find a guaranteed bracket."""
    return _pure_Tsat(P, 'ammonia'), _pure_Tsat(P, 'water')


# ---- optional Patek & Klomfar (1995) explicit correlations -----------------
_PK_B = [(0, 0, 0.322302e1), (0, 1, -0.384206), (0, 2, 0.460965e-1),
         (0, 3, -0.378945e-2), (0, 4, 0.135610e-3), (1, 0, 0.487755),
         (1, 1, -0.120108), (1, 2, 0.106154e-1), (2, 3, -0.533589e-3),
         (4, 0, 0.785041e1), (5, 0, -0.115941e2), (5, 1, -0.523150e-1),
         (6, 0, 0.489596e1), (13, 1, 0.421059e-1)]
_PK_D = [(0, 0, 0.324004e1), (0, 1, -0.395920), (0, 2, 0.435624e-1),
         (0, 3, -0.218943e-2), (1, 0, -0.143526e1), (1, 1, 0.105256e1),
         (1, 2, -0.719281e-1), (2, 0, 0.122362e2), (2, 1, -0.224368e1),
         (3, 0, -0.201780e2), (3, 1, 0.110834e1), (4, 0, 0.145399e2),
         (4, 2, 0.644312), (5, 0, -0.221246e1), (5, 2, -0.756266),
         (6, 0, -0.135529e1), (7, 2, 0.183541)]
_PK_T0, _PK_P0 = 100.0, 20.0     # K, bar (2 MPa)


def _pk_bubble_T(P, z):
    x = to_mole(z)
    L = np.log(_PK_P0 / P)
    return _PK_T0 * sum(a * (1 - x) ** m * L ** n for m, n, a in _PK_B)


def _pk_dew_T(P, z):
    y = to_mole(z)
    L = np.log(_PK_P0 / P)
    return _PK_T0 * sum(a * (1 - y) ** (m / 4) * L ** n for m, n, a in _PK_D)


# --------------------------------------------------------------------------
# state object
# --------------------------------------------------------------------------
@dataclass
class State:
    T: float; P: float; z: float; Q: float
    x: float; y: float
    h: float; s: float; v: float; u: float
    phase: str

    @property
    def rho(self):
        return 1.0 / self.v

    def as_dict(self):
        return asdict(self)

    def __repr__(self):
        return (f'State(T={self.T:.3f} K, P={self.P:.4g} bar, z={self.z:.4f}, '
                f'Q={self.Q:.4f}, x={self.x:.4f}, y={self.y:.4f}, '
                f'h={self.h:.3f} kJ/kg, s={self.s:.5f} kJ/kg/K, '
                f'v={self.v:.6g} m3/kg, {self.phase})')


def _assemble(T, P, z, warn=True):
    """TP-flash: the primitive every other flash reduces to."""
    check_range(T, P, warn)
    zx = to_mole(z)

    # pure components have no composition freedom -> classify by Tsat directly
    if zx <= _TINY or zx >= 1 - _TINY:
        fl = 'ammonia' if zx > 0.5 else 'water'
        try:
            Ts = _pure_Tsat(P, fl)
        except AmmoniaWaterError:
            Ts = None
        if Ts is None:
            ph = 'vapor'
        else:
            ph = 'liquid' if T < Ts else 'vapor'
        pr = props_single_phase(T, P, z, ph)
        return State(T, P, z, 0.0 if ph == 'liquid' else 1.0, z, z,
                     pr['h'], pr['s'], pr['v'], pr['u'], ph)

    # outside the pure-component saturation temperatures no mixture can split
    Ta, Tw = _envelope(P)
    if T <= Ta or T >= Tw:
        ph = 'liquid' if T <= Ta else 'vapor'
        pr = props_single_phase(T, P, z, ph)
        return State(T, P, z, 0.0 if ph == 'liquid' else 1.0, z, z,
                     pr['h'], pr['s'], pr['v'], pr['u'], ph)

    xl = _x_at_TP(T, P)
    yv = _y_from_x(T, P, xl)
    Xl, Yv = to_mass(xl), to_mass(yv)

    if not (Xl - 1e-9 <= z <= Yv + 1e-9):
        ph = 'liquid' if z < Xl else 'vapor'
        pr = props_single_phase(T, P, z, ph)
        return State(T, P, z, 0.0 if ph == 'liquid' else 1.0, z, z,
                     pr['h'], pr['s'], pr['v'], pr['u'], ph)

    Q = 0.0 if abs(Yv - Xl) < 1e-12 else (z - Xl) / (Yv - Xl)
    Q = min(max(Q, 0.0), 1.0)
    pl = props_single_phase(T, P, Xl, 'liquid')
    pg = props_single_phase(T, P, Yv, 'vapor')
    h = (1 - Q) * pl['h'] + Q * pg['h']
    s = (1 - Q) * pl['s'] + Q * pg['s']
    v = (1 - Q) * pl['v'] + Q * pg['v']
    return State(T, P, z, Q, Xl, Yv, h, s, v, h - P * v * 100.0, 'two-phase')


def _pure_two_phase(T, P, z, Q):
    """Saturated state of a pure component, where T and P are not independent
    and the quality must be supplied explicitly."""
    pl = props_single_phase(T, P, z, 'liquid')
    pg = props_single_phase(T, P, z, 'vapor')
    h = (1 - Q) * pl['h'] + Q * pg['h']
    s_ = (1 - Q) * pl['s'] + Q * pg['s']
    v = (1 - Q) * pl['v'] + Q * pg['v']
    ph = 'liquid' if Q <= 0 else ('vapor' if Q >= 1 else 'two-phase')
    return State(T, P, z, Q, z, z, h, s_, v, h - P * v * 100.0, ph)


def _is_pure(z):
    x = to_mole(z)
    return x <= _TINY or x >= 1 - _TINY


def _flash_TQ(T, z, Q, warn):
    """Find the pressure at which quality Q occurs at temperature T."""
    if _is_pure(z):
        P = _solve_bracketed(
            lambda P: _pure_Tsat(P, 'ammonia' if to_mole(z) > 0.5 else 'water') - T,
            0.02, 200.0, n=60, name='saturation pressure', log=True)
        check_range(T, P, warn)
        return _pure_two_phase(T, P, z, min(max(Q, 0.0), 1.0))
    # bubble/dew pressures bracket the two-phase region
    Pb = bubble_P(T, z)
    Pd = dew_P(T, z)
    lo, hi = min(Pb, Pd), max(Pb, Pd)
    if Q <= 1e-9:
        P = Pb
    elif Q >= 1 - 1e-9:
        P = Pd
    else:
        P = _solve_bracketed(lambda P: _assemble(T, P, z, warn=False).Q - Q,
                             lo * (1 + 1e-9), hi * (1 - 1e-9), n=60,
                             name='pressure', log=True)
    return _assemble(T, P, z, warn)


def _flash_PQ(P, z, Q, warn):
    if _is_pure(z):
        T = _pure_Tsat(P, 'ammonia' if to_mole(z) > 0.5 else 'water')
        check_range(T, P, warn)
        return _pure_two_phase(T, P, z, min(max(Q, 0.0), 1.0))
    Tb = bubble_T(P, z)
    Td = dew_T(P, z)
    if Q <= 1e-9:
        T = Tb
    elif Q >= 1 - 1e-9:
        T = Td
    else:
        T = _solve_bracketed(lambda T: _assemble(T, P, z, warn=False).Q - Q,
                             min(Tb, Td) + 1e-7, max(Tb, Td) - 1e-7,
                             n=40, name='temperature')
    return _assemble(T, P, z, warn)


_REL_TOL = {'h': 1e-4, 's': 1e-5, 'v': 1e-6, 'u': 1e-4}


def _pick_root(roots, span, key, target, unknown, warn):
    """Resolve a 1-D flash result, complaining when the inputs are ambiguous."""
    scale = max(abs(target), _REL_TOL.get(key, 1e-4))
    if not roots:
        if span < scale * 1e-3:
            raise AmmoniaWaterError(
                f'the given inputs do not determine the state: {key} varies by only '
                f'{span:.3g} over the whole {unknown} range, so {key} is effectively '
                f'independent of {unknown} here. Use a different input pair '
                f'(P,T or P,Q are always well posed).')
        raise AmmoniaWaterError(f'no solution for {key}={target:g}; '
                                f'attainable range spans {span:.4g}')
    if len(roots) > 1 and warn:
        warnings.warn(
            f'{len(roots)} states satisfy these inputs ({unknown} = '
            + ', '.join(f'{r:.4g}' for r in roots)
            + f'); returning the first. The pair (given, {key}) is not unique in '
              f'this region.', RuntimeWarning, stacklevel=4)
    return roots[0]


def _flash_Pk(P, z, key, target, warn):
    """(P, h) or (P, s) or (P, v) flash -- outer 1-D solve on temperature."""
    if _is_pure(z):
        st = _pure_saturated_pk(P, z, key, target, warn)
        if st is not None:
            return st

    def resid(T):
        return getattr(_assemble(T, P, z, warn=False), key) - target

    Ta, Tw = _envelope(P)
    roots, span = _all_roots(resid, max(205.0, Ta - 130.0), Tw + 220.0, n=150)
    T = _pick_root(roots, span, key, target, 'temperature', warn)
    return _assemble(T, P, z, warn)


def _pure_saturated_pk(P, z, key, target, warn):
    """For a pure component, check whether (P, key) lands inside the two-phase
    dome; if so the quality follows from the lever rule."""
    if key not in ('h', 's', 'v', 'u'):
        return None
    T = _pure_Tsat(P, 'ammonia' if to_mole(z) > 0.5 else 'water')
    lo = _pure_two_phase(T, P, z, 0.0)
    hi = _pure_two_phase(T, P, z, 1.0)
    a, b = getattr(lo, key), getattr(hi, key)
    if not (min(a, b) <= target <= max(a, b)):
        return None
    Q = 0.0 if abs(b - a) < 1e-14 else (target - a) / (b - a)
    check_range(T, P, warn)
    return _pure_two_phase(T, P, z, min(max(Q, 0.0), 1.0))


def _flash_Tk(T, z, key, target, warn):
    """(T, h) or (T, s) or (T, v) flash -- outer 1-D solve on pressure."""
    def resid(P):
        return getattr(_assemble(T, P, z, warn=False), key) - target

    roots, span = _all_roots(resid, 0.02, 190.0, n=150, log=True)
    P = _pick_root(roots, span, key, target, 'pressure', warn)
    return _assemble(T, P, z, warn)


_ALIAS = {'T': 'T', 'TEMPERATURE': 'T',
          'P': 'P', 'PRESSURE': 'P',
          'H': 'h', 'HMASS': 'h', 'ENTHALPY': 'h',
          'S': 's', 'SMASS': 's', 'ENTROPY': 's',
          'U': 'u', 'UMASS': 'u',
          'V': 'v', 'VMASS': 'v', 'VOLUME': 'v',
          'D': 'rho', 'DMASS': 'rho', 'RHO': 'rho', 'DENSITY': 'rho',
          'Q': 'Q', 'QUALITY': 'Q',
          'X': 'x', 'Y': 'y', 'Z': 'z', 'PHASE': 'phase'}


def state(z, warn=True, **kwargs):
    """Resolve a full state from the NH3 mass fraction z and any valid pair.

    Accepted keys: T, P, H, S, U, V, D, Q  (case-insensitive).

        state(z=0.4, T=350, P=10)
        state(z=0.4, P=10, Q=0)        # bubble point
        state(z=0.4, P=10, H=500)
        state(z=0.4, T=350, S=1.2)
    """
    if not (0.0 <= z <= 1.0):
        raise AmmoniaWaterError(f'z must lie in [0, 1], got {z}')
    given = {}
    for k, val in kwargs.items():
        kk = _ALIAS.get(k.upper())
        if kk is None:
            raise AmmoniaWaterError(f'unknown input {k!r}')
        given[kk] = float(val)
    if len(given) != 2:
        raise AmmoniaWaterError(f'need exactly 2 inputs besides z, got {len(given)}')

    if 'rho' in given:                       # density -> specific volume
        given['v'] = 1.0 / given.pop('rho')
    keys = set(given)

    if keys == {'T', 'P'}:
        return _assemble(given['T'], given['P'], z, warn)
    if keys == {'P', 'Q'}:
        return _flash_PQ(given['P'], z, given['Q'], warn)
    if keys == {'T', 'Q'}:
        return _flash_TQ(given['T'], z, given['Q'], warn)
    if 'P' in keys:
        other = (keys - {'P'}).pop()
        return _flash_Pk(given['P'], z, other, given[other], warn)
    if 'T' in keys:
        other = (keys - {'T'}).pop()
        return _flash_Tk(given['T'], z, other, given[other], warn)
    raise AmmoniaWaterError(f'input pair {sorted(keys)} is not supported; '
                            f'include T or P')


def PropsSI(output, n1, v1, n2, v2, z, warn=True):
    """CoolProp-style single-value accessor.

        PropsSI('T', 'P', 10, 'Q', 0, 0.4)   ->  bubble temperature [K]
        PropsSI('H', 'T', 350, 'P', 10, 0.4) ->  enthalpy [kJ/kg]
    """
    st = state(z, warn=warn, **{n1: v1, n2: v2})
    key = _ALIAS.get(output.upper())
    if key is None:
        raise AmmoniaWaterError(f'unknown output {output!r}')
    return getattr(st, key)


# --------------------------------------------------------------------------
# ThermoSim integration: a "fluid" identifier that carries composition
# --------------------------------------------------------------------------
class NH3H2OFluid:
    """Identifies an NH3-H2O stream of a fixed overall NH3 mass fraction.

    Pass this wherever ThermoSim expects a CoolProp fluid-name string, e.g.
    ``Model.add_point(NH3H2OFluid(0.4), 'S1', P=..., T=...)``. `Prop`
    dispatches on `isinstance(fluid, NH3H2OFluid)` (mirroring how it already
    dispatches on the Therminol-66 fluid-name aliases) and routes property
    calculations through this module instead of CoolProp.

    Composition is bundled into the fluid identifier -- not passed as a
    separate argument -- specifically so that every existing component which
    does ``self.Model.Prop(self.In.fluid, P=..., S=...)`` (Pump, Turbine,
    Compressor, Expansion_valve, and the HeatExchanger helpers) keeps working
    unchanged: whatever `z` the inlet point carries is forwarded automatically
    because it lives on the object those call sites already pass through.

    Only components that deliberately change composition (a Separator
    splitting a two-phase mixture into its equilibrium liquid/vapor
    compositions) construct a *new* `NH3H2OFluid` for their outlet.
    """

    __slots__ = ('z',)
    name = 'NH3H2O'

    def __init__(self, z):
        if not (0.0 <= z <= 1.0):
            raise ValueError(f"NH3H2OFluid: z must lie in [0, 1], got {z}")
        self.z = float(z)

    def __eq__(self, other):
        return isinstance(other, NH3H2OFluid) and self.z == other.z

    def __hash__(self):
        return hash(('NH3H2O', self.z))

    def __repr__(self):
        return f'NH3H2O(z={self.z:.4f})'

    __str__ = __repr__
