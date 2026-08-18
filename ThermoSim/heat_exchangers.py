"""
heat_exchangers.py
------------------
Generalised HeatExchanger solver supporting:
  - 0, 1, 2, and 3 unknown variables from {H_hi, H_ho, H_ci, H_co, mh, mc}
  - All solvable combinations via an anchor + brentq engine
  - HEX types: double_pipe, Condenser, Evaporator, SimpleHEX
  - TES (Thermal Energy Storage)

Solvability
-----------
The authoritative reference is ``docs/HEX_Logic_Table.xlsx``; the counts below
match its ✅/❌ columns, the CHANGELOG and ``_SOLVABLE_2_UNKNOWN`` below.

  2-unknown: 10 of C(6,2) = 15 combinations are solvable.
             Unsolvable: cases 1, 2 (same-side enthalpy pair → Q undefined),
             8, 9 (knowns confined to one end → unreliable) and
             15 (mh/mc ratio fixed → infinitely many solutions).

  3-unknown:  8 of C(6,3) = 20 combinations are solvable.
             Strategy: anchor one unknown via the PPT diagonal, then hand the
             remaining pair to the 2-unknown solver.  Two pinch points are
             found in total, so the T-profile scan excludes the anchored
             index.  The other 12 either have no valid anchor (rows 11-14) or
             reduce to an unsolvable 2-unknown case (rows 1, 2, 5, 6 → cases
             9/8; rows 17-20 → case 15).

Earlier releases documented "10/15" and "18/20" here; the 3-unknown figure was
never correct.
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import CoolProp
import CoolProp.CoolProp as CP
from scipy.optimize import brentq, differential_evolution

from .base_component import Component


# --------------------------------------------------------------------------- #
#  Lower temperature bound for root-finding brackets
# --------------------------------------------------------------------------- #
_AS_CACHE = {}


def _safe_T_min(fluid, P, margin=0.05, max_steps=40):
    """
    Lowest temperature CoolProp will actually evaluate for `fluid` at `P`.

    The pinch solver brackets its search at the coldest valid state, and the
    obvious choice -- ``PropsSI('TMIN', fluid)`` -- is not always valid.
    TMIN is a single number, but the melting line Tmelt(p) varies with
    pressure and can rise ABOVE it, at which point CoolProp refuses the
    evaluation:

        For now, we don't support T [59.75 K] below Tmelt(p) [59.9275 K]

    That made every 2- and 3-unknown pinch problem fail for Air, Nitrogen,
    CO2 and Oxygen -- i.e. for every gas-turbine, HRSG and recuperator model.

    Water is the opposite case: its melting line FALLS with pressure, because
    ice is less dense than water. Taking the max of TMIN and Tmelt(p)
    therefore handles both directions with one expression.

    Not every fluid has a melting line. Refrigerants have none fitted, and
    CO2 below its triple-point pressure sublimes rather than melts. Both
    raise, and both correctly fall back to TMIN.

    The bound is finally confirmed by evaluation and walked upward if that
    fails, so this never returns a temperature CoolProp will reject.
    """
    T_lo = CP.PropsSI('TMIN', fluid)
    try:
        if fluid not in _AS_CACHE:
            _AS_CACHE[fluid] = CoolProp.AbstractState('HEOS', fluid)
        T_melt = _AS_CACHE[fluid].melting_line(CoolProp.iT, CoolProp.iP, P)
        T_lo = max(T_lo, T_melt)
    except Exception:
        pass                       # no melting line available; TMIN it is

    T = T_lo + margin
    for _ in range(max_steps):
        try:
            CP.PropsSI('H', 'T', T, 'P', P, fluid)
            return T
        except Exception:
            T += max(margin, 0.01 * T_lo)
    raise ValueError(
        f"_safe_T_min: no valid lower temperature bound for {fluid} "
        f"at P = {P:.0f} Pa"
    )



# ─────────────────────────────────────────────────────────────────────────────
#  Internal helpers
# ─────────────────────────────────────────────────────────────────────────────
#: Value returned by a residual function when the trial state is
#: thermodynamically invalid.  Using an explicit sentinel (rather than a
#: signed magnitude such as ``1e10 * sign_factor``) means the bracket
#: inspection in ``_fun_sign_fac`` can recognise a penalty reliably; with a
#: signed magnitude, a penalty of ``-1e10`` never matched the ``== 1e10``
#: tests and the sign heuristic silently stopped working after its first use.
_PENALTY = 1e10


def robust_solver(func, a, b, tol=1e-4):
    """
    Fallback root-finder used when ``brentq`` cannot bracket a root.

    Minimises ``|func(x)|`` with Differential Evolution over ``[a, b]``.
    Slower than brentq but does not require a sign change at the ends, which
    matters here because the residual is a discontinuous min(ΔT) surface.

    Raises
    ------
    ValueError
        If the best point found is not close enough to a genuine root.
    """
    if b <= a:
        raise ValueError(
            f"robust_solver: empty search bracket [{a:.6g}, {b:.6g}]."
        )

    def _objective(x):
        """
        Scalarising wrapper around the residual.

        differential_evolution passes the trial point as an ndarray of
        shape (1,), but the residual functions expect a plain float: handing
        them an array makes CoolProp raise, the bare-except path returns the
        scalar _PENALTY, and the mix of scalar and array returns makes SciPy
        abort with "The map-like callable must be of the form f(func,
        iterable)".  Unwrapping here keeps every residual scalar-only.
        """
        val = func(float(np.ravel(x)[0]))
        return abs(float(np.ravel(val)[0]))

    res = differential_evolution(
        _objective,
        bounds=[(a, b)],
        tol=tol,
        polish=True,      # quick local search at the end for precision
        popsize=10,       # increase if the function is extremely jagged
    )

    if res.fun > 1e-2:
        raise ValueError(
            f"robust_solver: no root found in [{a:.6g}, {b:.6g}]. "
            f"Best residual = {res.fun:.6g} K at x = {res.x[0]:.6g}."
        )

    return res.x[0]


def _h_from_T(model, fluid, state_name, T, P, mflow=None):
    """Return a Prop state computed from (T, P)."""
    return model.Prop(fluid, StatePointName=state_name, T=T, P=P,
                      Mass_flowrate=mflow)


def _h_from_H(model, fluid, state_name, H, P, mflow=None):
    """Return a Prop state computed from (H, P)."""
    return model.Prop(fluid, StatePointName=state_name, H=H, P=P,
                      Mass_flowrate=mflow)


def _T_sat(model, fluid, P, Q_val):
    """Saturation temperature at given pressure and vapour quality."""
    return model.Prop(fluid, StatePointName='_sat', P=P, Q=Q_val).T


# ─────────────────────────────────────────────────────────────────────────────
#  HeatExchanger
# ─────────────────────────────────────────────────────────────────────────────
class HeatExchanger(Component):
    """
    Generalised counter-flow heat exchanger.

    Parameters
    ----------
    Model          : simulation model object (holds Point dict and Prop factory)
    ID             : string identifier
    PPT            : pinch-point temperature difference [K]
    HEX_type       : 'double_pipe' | 'Condenser' | 'Evaporator' | 'SimpleHEX'
    HeatAdded      : True (heat source) | False (heat sink) | None (internal HEX)
    Hot_In_state   : key into Model.Point for hot inlet  (None for SimpleHEX one-side)
    Hot_Out_state  : key into Model.Point for hot outlet
    Cold_In_state  : key into Model.Point for cold inlet
    Cold_Out_state : key into Model.Point for cold outlet
    UA             : optional overall conductance [W/K]
    effectiveness  : optional NTU-ε value (0-1)
    Q              : optional heat duty [W]
    div_N          : number of segments for pinch profile (default 200)
    PPT_graph      : if True, display pinch temperature profile plot
    Calculate      : if True, call Cal() immediately on construction
    """

    #: The 10 solvable 2-unknown combinations, matching the ✅ column of
    #: HEX_Logic_Table.xlsx (sheet "2-Unknown Cases").  The five omitted
    #: combinations are genuinely unsolvable, not merely unimplemented:
    #:   {H_hi, H_ho}  case 1  — both hot enthalpies unknown, Q undefined
    #:   {H_ci, H_co}  case 2  — both cold enthalpies unknown, Q undefined
    #:   {H_hi, H_co}  case 8  — knowns sit at one end only; unreliable
    #:   {H_ho, H_ci}  case 9  — knowns sit at one end only; unreliable
    #:   {mh,   mc}    case 15 — mh/mc ratio is fixed → infinitely many roots
    _SOLVABLE_2_UNKNOWN = frozenset({
        frozenset({'H_hi', 'mh'}),    # case 3
        frozenset({'H_ho', 'mh'}),    # case 4
        frozenset({'H_ci', 'mc'}),    # case 5
        frozenset({'H_co', 'mc'}),    # case 6
        frozenset({'H_hi', 'H_ci'}),  # case 7
        frozenset({'H_ho', 'H_co'}),  # case 10
        frozenset({'H_hi', 'mc'}),    # case 11
        frozenset({'H_ho', 'mc'}),    # case 12
        frozenset({'H_ci', 'mh'}),    # case 13
        frozenset({'H_co', 'mh'}),    # case 14
    })

    def __init__(self, Model, ID, PPT, HEX_type, HeatAdded,
                 Hot_In_state, Hot_Out_state,
                 Cold_In_state, Cold_Out_state,
                 UA=None, effectiveness=None, Q=None,
                 div_N=200, PPT_graph=False, Calculate=False,print_residual = False):

        self.HeatAdded      = HeatAdded
        self.Hot_In_state   = Hot_In_state
        self.Hot_Out_state  = Hot_Out_state
        self.Cold_In_state  = Cold_In_state
        self.Cold_Out_state = Cold_Out_state
        self.HEX_type       = HEX_type
        self.div_N          = div_N
        self.PPT            = PPT
        self.PPT_graph      = PPT_graph
        self.Hot_to_Cold    = None
        self.LMTD           = None
        self.Q              = Q
        self.UA             = UA
        self.effectiveness  = effectiveness
        self.Hot_Mass_flowrate  = None
        self.Cold_Mass_flowrate = None
        self.print_residual = print_residual

        if HeatAdded not in (True, False, None):
            raise ValueError(
                "HeatAdded must be True, False, or None."
            )

        # Registration, Solution_Status and Ex_D are handled by the base
        # class; v3.2.1 duplicated all of it here and never called super().
        super().__init__(Model, ID, Calculate)

    # ── helpers ──────────────────────────────────────────────────────────────

    @staticmethod
    def _is_two_phase(point):
        try:
            q = CP.PropsSI('Q', 'H', point.H, 'P', point.P, point.fluid)
            return 0.0 < q < 1.0
        except Exception:
            return False

    def _min_H(self, pt, margin=1e-3, T_floor=None):
        """
        Lower enthalpy bound for a root-finding bracket on `pt`.

        Two limits apply, and the tighter one wins.

        The first is what CoolProp will evaluate at all -- see
        ``_safe_T_min``.  On its own that is a very weak bound for gases:
        air is valid down to about 60 K, several hundred kelvin below any
        state a heat exchanger will actually visit.  Bracketing over a range
        that wide leaves the pinch residual flat and sign-free across most
        of it, ``brentq`` cannot bracket a root, and the solver falls back to
        differential evolution -- correct, but around fifty times slower.

        The second is physical.  In counter-flow neither stream can cross
        the other, so a hot stream never cools below the cold inlet and a
        cold stream never leaves below its own inlet.  Passing that
        temperature as ``T_floor`` tightens the bracket to the region a
        solution can actually occupy, which is both faster and less prone to
        converging on a spurious root.
        """
        T = _safe_T_min(pt.fluid, pt.P)

        # Default floor: the cold inlet.  It is the coldest state anywhere in
        # a counter-flow exchanger, so no bracket needs to reach below it.
        if T_floor is None:
            CI = getattr(self, 'Cold_In', None)
            if CI is not None and getattr(CI, 'T', None) is not None:
                T_floor = CI.T

        if T_floor is not None and T_floor > T:
            T = T_floor
        return self.Model.Prop(pt.fluid, StatePointName='_lo_bound',
                               T=T, P=pt.P).H + margin

    def _resolve_mass_flow(self, pt_in, pt_out, side):
        m_in  = pt_in.Mass_flowrate  if pt_in  else None
        m_out = pt_out.Mass_flowrate if pt_out else None
        if m_in is None and m_out is None:
            return None
        if m_out is None:
            if pt_out: pt_out.Mass_flowrate = m_in
            return m_in
        if m_in is None:
            if pt_in: pt_in.Mass_flowrate = m_out
            return m_out
        if abs(m_in - m_out) / max(abs(m_in), 1e-12) < 1e-6:
            return m_in
        raise ValueError(
            f"[{self.ID}] {side}-side mass-flow mismatch: "
            f"inlet={m_in} vs outlet={m_out} kg/s."
        )

    def _check_energy_balance(self, tol=0.0001):
        if self.Hot_Mass_flowrate is None or self.Cold_Mass_flowrate is None:
            return
        Q_h = (self.Hot_In.H  - self.Hot_Out.H)  * self.Hot_Mass_flowrate
        Q_c = (self.Cold_Out.H - self.Cold_In.H) * self.Cold_Mass_flowrate
        ref = max(abs(Q_h), abs(Q_c), 1.0)
        if abs(Q_h - Q_c) / ref > tol:
            raise ValueError(
                f"[{self.ID}] Energy imbalance: Q_hot={Q_h:.2f} W, "
                f"Q_cold={Q_c:.2f} W (err={abs(Q_h-Q_c)/ref*100:.2f} %)."
            )
        if Q_c < 0 and Q_h < 0:
            raise ValueError(
                f"[{self.ID}] Heat flows the wrong way: Q_hot={Q_h:.2f} W and "
                f"Q_cold={Q_c:.2f} W are both negative, i.e. the 'cold' stream "
                f"is hotter than the 'hot' stream. Check which streams are "
                f"assigned to Hot_* and Cold_*."
            )

    def _writeback_point(self, state_key, point):
        if state_key:
            self.Model.Point[state_key] = point

    def _writeback_mflow_hot(self, mh):
        self.Hot_Mass_flowrate = mh
        self.Hot_In.Mass_flowrate  = mh
        self.Hot_Out.Mass_flowrate = mh
        self._writeback_point(self.Hot_In_state,  self.Hot_In)
        self._writeback_point(self.Hot_Out_state, self.Hot_Out)

    def _writeback_mflow_cold(self, mc):
        self.Cold_Mass_flowrate = mc
        self.Cold_In.Mass_flowrate  = mc
        self.Cold_Out.Mass_flowrate = mc
        self._writeback_point(self.Cold_In_state,  self.Cold_In)
        self._writeback_point(self.Cold_Out_state, self.Cold_Out)

    # ── pinch profile ─────────────────────────────────────────────────────────

    def _build_profile(self, Q_total, dP_h, dP_c,
                       H_h0, H_c0, P_h0, P_c0,
                       mh, mc,
                       direction='cold_to_hot',
                       exclude_idx=None):
        """
        Build discretised T-profile arrays.

        direction='cold_to_hot'  : marching from hot outlet / cold inlet
        direction='hot_to_cold'  : marching from hot inlet  / cold outlet

        exclude_idx : if set, the residual skips that index when computing min(ΔT).
                      Used in 3-unknown solver to avoid the anchored pinch point.

        Returns Th, Tc, dT arrays of length N+1.
        """
        N  = self.div_N
        q  = Q_total / N
        Th = np.zeros(N + 1)
        Tc = np.zeros(N + 1)
        dT = np.zeros(N + 1)

        Th[0] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h0',
                                H=H_h0, P=P_h0).T
        Tc[0] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c0',
                                H=H_c0, P=P_c0).T
        dT[0] = Th[0] - Tc[0]
        
        for n in range(1, N + 1):
            if direction == 'cold_to_hot':
                H_h = H_h0 + q * n / mh
                H_c = H_c0 + q * n / mc
                Ph  = P_h0 + dP_h * n / N
                Pc  = P_c0 - dP_c * n / N
            else:   # direction == 'hot_to_cold'
                H_h = H_h0 - q * n / mh
                H_c = H_c0 - q * n / mc
                Ph  = P_h0 - dP_h * n / N
                Pc  = P_c0 + dP_c * n / N
            
            Th[n] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h',
                                    H=H_h, P=Ph).T
            Tc[n] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c',
                                    H=H_c, P=Pc).T
            dT[n] = Th[n] - Tc[n]

        if exclude_idx is not None:
            mask = np.ones(N + 1, dtype=bool)
            if (exclude_idx == 0 and direction == 'cold_to_hot') or (exclude_idx != 0 and direction == 'hot_to_cold'):
                mask[0:int(N*0.1)] = False
            elif (exclude_idx != 0 and direction == 'cold_to_hot') or (exclude_idx == 0 and direction == 'hot_to_cold'):
                mask[int(N*0.9):] = False   # slice to the END: index N is the hot-inlet
                #                     anchor point and must be excluded too
            return Th, Tc, dT, np.min(dT[mask])

        return Th, Tc, dT, np.min(dT)

    def _residual(self, Q_total, dP_h, dP_c,
                  H_h0, H_c0, P_h0, P_c0,
                  mh, mc, direction='cold_to_hot',
                  exclude_idx=None):
        """Return min(ΔT) − PPT for use as a brentq residual."""
        _, _, _, min_dT = self._build_profile(
            Q_total, dP_h, dP_c,
            H_h0, H_c0, P_h0, P_c0,
            mh, mc, direction, exclude_idx)
        if self.print_residual:
            print(f'solving {self.ID}. Residual value: {min_dT - self.PPT}')
        return min_dT - self.PPT

    # ── pinch check helpers ───────────────────────────────────────────────────

    def _check_PPT_ends(self):
        """Verify PPT is not violated at either end after solving."""
        if self.Hot_In.T is not None and self.Cold_Out.T is not None:
            dT_hi = self.Hot_In.T  - self.Cold_Out.T
        else:
            dT_hi = None

        if self.Hot_Out.T is not None and self.Cold_In.T is not None:
            dT_ho = self.Hot_Out.T - self.Cold_In.T
        else:
            dT_ho = None

        # The pinch is located by a root-finder working on enthalpy with a
        # RELATIVE tolerance (brentq rtol=1e-4), so the temperature residual it
        # leaves behind scales with the size of the problem.  A fixed absolute
        # allowance therefore rejects perfectly good solutions on large pinches:
        # a 30 K PPT converging to 29.986 K is 0.05 % out, i.e. numerical noise,
        # yet a flat 1e-2 K tolerance called it a physical violation.  Scale the
        # allowance with PPT, keeping the old strictness for small pinches.
        # (Same reasoning as the Solve() convergence tolerance in 3.2.2.)
        tol = max(1e-2, 1e-3 * self.PPT)

        for name, val in [("hot-inlet/cold-outlet", dT_hi),
                          ("hot-outlet/cold-inlet", dT_ho)]:
            if val is not None:
                if val <= self.PPT - tol:
                    raise ValueError(
                        f"[{self.ID}] PPT violation at {name}: "
                        f"ΔT={val:.3f} K < PPT={self.PPT:.3f} K."
                    )

    # ─────────────────────────────────────────────────────────────────────────
    #  2-UNKNOWN SOLVER ENGINE
    # ─────────────────────────────────────────────────────────────────────────

    def _solve_2_unknown(self, unknowns, dP_h, dP_c, exclude_idx=None):
        """
        Generalised 2-unknown solver.

        unknowns : frozenset of 2 strings from
                   {'H_hi','H_ho','H_ci','H_co','mh','mc'}

        exclude_idx : passed to _build_profile to skip the anchor pinch point
                      when called from the 3-unknown solver.

        The solver always uses 'cold_to_hot' direction (marching from hot outlet
        / cold inlet) unless only hot-outlet/cold-inlet are defined.

        Strategy per case:
          Free variable (brentq) → derive dependent variable from energy balance
          → compute T-profile → return min(ΔT) − PPT.
        """
        HI = self.Hot_In;   HO = self.Hot_Out
        CI = self.Cold_In;  CO = self.Cold_Out
        mh = self.Hot_Mass_flowrate
        mc = self.Cold_Mass_flowrate
        PPT = self.PPT
        N   = self.div_N
        mdl = self.Model
        
        # ── Unsolvable guards (HEX_Logic_Table.xlsx, "2-Unknown Cases") ─────
        _unsolvable = {
            frozenset({'H_hi', 'H_ho'}):
                "both hot-side enthalpies are unknown, so Q is undefined on "
                "that side (case 1)",
            frozenset({'H_ci', 'H_co'}):
                "both cold-side enthalpies are unknown, so Q is undefined on "
                "that side (case 2)",
            frozenset({'H_hi', 'H_co'}):
                "the known enthalpies both sit at the hot-outlet/cold-inlet "
                "end, so a solution may exist but is not reliable (case 8)",
            frozenset({'H_ho', 'H_ci'}):
                "the known enthalpies both sit at the hot-inlet/cold-outlet "
                "end, so a solution may exist but is not reliable (case 9)",
            frozenset({'mh', 'mc'}):
                "mh and mc only ever appear as a fixed ratio, so there are "
                "infinitely many solutions (case 15)",
        }
        if unknowns in _unsolvable:
            raise ValueError(
                f"[{self.ID}] Unsolvable 2-unknown combination "
                f"{set(unknowns)}: {_unsolvable[unknowns]}."
            )
        if unknowns not in self._SOLVABLE_2_UNKNOWN:
            raise ValueError(
                f"[{self.ID}] Unrecognised 2-unknown combination "
                f"{set(unknowns)}. Solvable combinations are: "
                f"{sorted(set(u) for u in self._SOLVABLE_2_UNKNOWN)}."
            )

        # ── helper: profile residual with cold_to_hot direction ──────────────
        def _res_cth(Q_sp, mh_, mc_):
            HO = self.Hot_Out
            CI = self.Cold_In

            return self._residual(
                Q_sp, dP_h, dP_c,
                HO.H, CI.H, HO.P, CI.P,
                mh_, mc_, 'cold_to_hot', exclude_idx)

        def _res_htc(Q_sp, mh_, mc_):
            HI = self.Hot_In
            CO = self.Cold_Out
            return self._residual(
                Q_sp, dP_h, dP_c,
                HI.H, CO.H, HI.P, CO.P,
                mh_, mc_, 'hot_to_cold', exclude_idx)
        
        def _fun_sign_fac(f, H_lo, H_hi_bound):
            """
            Choose the sign of the invalid-state penalty so that the penalty
            sits on the OPPOSITE side of zero from the valid end, giving
            brentq a sign change to bracket.

            The penalty is compared against the bare ``_PENALTY`` magnitude.
            v3.2.1 compared against ``1e10`` while ``f`` returned
            ``1e10 * sign_factor``, so once the sign flipped to -1 none of
            the tests could ever match again and the heuristic became inert.
            """
            nonlocal sign_factor
            sign_factor = 1                      # probe with a positive penalty
            Hi = f(H_hi_bound)
            Lo = f(H_lo)

            lo_bad = (Lo == _PENALTY)
            hi_bad = (Hi == _PENALTY)

            if lo_bad and not hi_bad:
                return -1 if Hi > 1e-7 else 1
            if hi_bad and not lo_bad:
                return -1 if Lo > 1e-7 else 1
            return 1


        # ════════════════════════════════════════════════════════════════════
        # CASES involving one mass-flow unknown + one enthalpy
        # ════════════════════════════════════════════════════════════════════
        sign_factor = 1
        # H_hi unknown, mh known, mc known  →  free: mh, derive H_hi
        if unknowns == frozenset({'H_hi', 'mh'}):
            Q_c = mc * (CO.H - CI.H)
            def f(H_hi):
                try:
                    mh_ = Q_c / (H_hi - HO.H)   
                    return _res_cth(Q_c, mh_, mc)
                except Exception:
                    return _PENALTY * sign_factor  # large penalty if state point is invalid
            H_lo = self._min_H(HI)
            H_hi_bound = self.Model.Prop(HI.fluid, StatePointName='_hi_bound', T=CP.PropsSI('TMAX', HI.fluid), P=HI.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:
                try:
                    H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    _H_lo = HO.H + 1e-3
                    sign_factor = _fun_sign_fac(f,_H_lo,H_hi_bound)
                    H_sol = brentq(f, _H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_sol={H_sol:.2f} J/kg")
            mh_sol = Q_c / (H_sol - HO.H)
            self.Hot_In = _h_from_H(mdl, HI.fluid, HI.StatePointName, H_sol, HI.P, mh_sol)
            self._writeback_mflow_hot(mh_sol)
            return

        if unknowns == frozenset({'H_ho', 'mh'}):
            Q_c = mc * (CO.H - CI.H)
            def f(H_ho):
                try:
                    self.Hot_Out = _h_from_H(mdl, HO.fluid, HO.StatePointName, H_ho, HO.P)
                    mh_ = Q_c / (HI.H - H_ho)
                    return _res_htc(Q_c, mh_, mc)
                except Exception:
                    return _PENALTY * sign_factor
            H_lo = self._min_H(HO, T_floor=CI.T if CI is not None else None)
            H_hi_bound = self.Model.Prop(HO.fluid, StatePointName='_ho_hi', T=CP.PropsSI('TMAX', HO.fluid), P=HO.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:
                try:
                    H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    _H_hi_bound = HI.H-1e-3
                    sign_factor = _fun_sign_fac(f,H_lo,_H_hi_bound)
                    H_sol = brentq(f, H_lo, _H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_sol={H_sol:.2f} J/kg")
            
            mh_sol = Q_c / (HI.H - H_sol)   
            self.Hot_Out = _h_from_H(mdl, HO.fluid, HO.StatePointName, H_sol, HO.P, mh_sol)
            self._writeback_mflow_hot(mh_sol)
            return

        if unknowns == frozenset({'H_ci', 'mc'}):
            Q_h = mh * (HI.H - HO.H)
            def f(H_ci):
                try:
                    self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H_ci, CI.P)
                    mc = Q_h / (CO.H - self.Cold_In.H)
                    return _res_htc(Q_h, mh, mc)
                except Exception:
                    return _PENALTY * sign_factor  # large penalty if state point is invalid
            H_lo = self._min_H(CI)
            H_hi_bound = self.Model.Prop(CI.fluid, StatePointName='_ci_hi', T=CP.PropsSI('TMAX', CI.fluid), P=CI.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:
                try:
                    H_ci = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    _H_hi_bound = CO.H-1e-3
                    sign_factor = _fun_sign_fac(f,H_lo,_H_hi_bound)
                    H_ci = brentq(f, H_lo, _H_hi_bound, xtol=1e-4, rtol=1e-4)

            except Exception:
                H_ci = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_ci={H_sol:.2f} J/kg")

            mc = Q_h/(CO.H - H_ci)
            self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H_ci, CI.P, mc)
            self._writeback_mflow_cold(mc)
            return

        if unknowns == frozenset({'H_co', 'mc'}):
            Q_h = mh * (HI.H - HO.H)
            def f(H_co):
                try:
                    self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H_co, CO.P)
                    mc_ = Q_h / (H_co - CI.H)
                    return _res_cth(Q_h, mh, mc_)
                except Exception:
                    return _PENALTY * sign_factor  # large penalty if state point is invalid
            H_lo = self._min_H(CO, T_floor=CI.T if CI is not None else None)
            H_hi_bound = self.Model.Prop(CO.fluid, StatePointName='_co_hi', T=CP.PropsSI('TMAX', CO.fluid), P=CO.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:
                try:
                    H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    # Narrow the upper bound to the cold-side enthalpy that
                    # corresponds to the hot inlet temperature.  Bounding a
                    # COLD enthalpy by the HOT enthalpy HI.H (as v3.2.1 did)
                    # is meaningless whenever the two streams are different
                    # fluids, because their enthalpy references differ.
                    _H_hi_bound = self.Model.Prop(
                        CO.fluid, StatePointName='_co_hi2',
                        P=CO.P, T=HI.T).H - 1e-3
                    sign_factor = _fun_sign_fac(f, H_lo, _H_hi_bound)
                    # v3.2.1 passed H_hi_bound here, i.e. it retried with the
                    # exact bracket that had just failed.
                    H_sol = brentq(f, H_lo, _H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_co={H_sol:.2f} J/kg")

            mc = Q_h / (H_sol - CI.H)
            self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H_sol, CO.P, mc)
            self._writeback_mflow_cold(mc)
            return

        if unknowns == frozenset({'H_hi', 'mc'}):
            dH_c = CO.H - CI.H
            def f(H_hi):
                try:
                    self.Hot_In = _h_from_H(mdl, HI.fluid, HI.StatePointName, H_hi, HI.P, mh)
                    Q_h = mh * (self.Hot_In.H - HO.H)
                    mc_ = Q_h / dH_c
                except Exception:
                    return _PENALTY * sign_factor  # large penalty if state point is invalid
                return _res_cth(Q_h, mh, mc_)
            H_lo = self._min_H(HI)
            H_hi_bound = self.Model.Prop(HI.fluid, StatePointName='_hi_bound', T=CP.PropsSI('TMAX', HI.fluid) , P=HI.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:
                try:
                    H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    _H_lo = HO.H + 1e-3
                    sign_factor = _fun_sign_fac(f,_H_lo,H_hi_bound)
                    H_sol = brentq(f, _H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_hi={H_sol:.2f} J/kg")
            
            self.Hot_In = _h_from_H(mdl, HI.fluid, HI.StatePointName, H_sol, HI.P, mh)
            Q_h  = mh * (self.Hot_In.H - HO.H)
            mc_sol = Q_h / dH_c
            self._writeback_mflow_cold(mc_sol)
            return

        if unknowns == frozenset({'H_ho', 'mc'}):
            dH_c = CO.H - CI.H
            def f(H_ho):
                try:
                    # self.Hot_Out = _h_from_T(mdl, HO.fluid, HO.StatePointName, T_ho, HO.P)
                    Q_h = mh * (HI.H - H_ho)
                    mc_ = Q_h / dH_c
                    return _res_htc(Q_h, mh, mc_)
                except Exception:
                    return _PENALTY * sign_factor  # large penalty if state point is invalid
            H_lo = self._min_H(HO, margin=1.0,
                               T_floor=CI.T if CI is not None else None)
            H_hi_bound = self.Model.Prop(HO.fluid, StatePointName='_ho_hi', T= CP.PropsSI('TMAX', HO.fluid), P=HO.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:
                try:
                    H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    # H_ho and HI.H are the same (hot) stream, so this bound
                    # is dimensionally sound: the outlet cannot exceed the
                    # inlet enthalpy for a cooling stream.
                    _H_hi_bound = HI.H - 1e-3
                    sign_factor = _fun_sign_fac(f, H_lo, _H_hi_bound)
                    H_sol = brentq(f, H_lo, _H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_ho={H_sol:.2f} J/kg")

            self.Hot_Out = _h_from_H(mdl, HO.fluid, HO.StatePointName, H_sol, HO.P,mh)
            Q_h  = mh * (HI.H - self.Hot_Out.H)
            mc_sol = Q_h / dH_c
            self._writeback_mflow_cold(mc_sol)
            return

        if unknowns == frozenset({'H_ci', 'mh'}):
            dH_h = HI.H - HO.H
            def f(H_ci):
                try:
                    Q_c = mc * (CO.H - H_ci)
                    mh_ = Q_c / dH_h
                    return _res_htc(Q_c, mh_, mc)
                except Exception:
                    return _PENALTY * sign_factor
            H_lo = self._min_H(CI)
            H_hi_bound = self.Model.Prop(CI.fluid, StatePointName='_lo_bound', T=CP.PropsSI('TMAX', CI.fluid), P=CI.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:                
                try:
                    H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    _H_hi_bound = CO.H-1e-3
                    sign_factor = _fun_sign_fac(f,H_lo,_H_hi_bound)
                    H_sol = brentq(f, H_lo, _H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:                
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_ci={H_sol:.2f} J/kg")
        
            Q_c  = mc * (CO.H - H_sol)
            mh_sol = Q_c / dH_h
            self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H_sol, CI.P,mc)
            self._writeback_mflow_hot(mh_sol)
            return

        if unknowns == frozenset({'H_co', 'mh'}):
            dH_h = HI.H - HO.H
            def f(H_co):
                try:
                    Q_c = mc * (H_co - CI.H)
                    mh_ = Q_c / dH_h
                    return _res_cth(Q_c, mh_, mc)
                except Exception:
                    return _PENALTY * sign_factor
            H_lo = self._min_H(CO, T_floor=CI.T if CI is not None else None)
            H_hi_bound = self.Model.Prop(CO.fluid, StatePointName='_hi_bound', T=CP.PropsSI('TMAX', CO.fluid), P=CO.P).H-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:      
                try:          
                    H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
                except Exception:
                    _H_lo = CI.H - 1e-3
                    sign_factor = _fun_sign_fac(f,_H_lo,H_hi_bound)
                    H_sol = brentq(f, _H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:                
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_ci={H_sol:.2f} J/kg")

            self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H_sol, CO.P,mc)
            Q_c  = mc * (self.Cold_Out.H - CI.H)
            mh_sol = Q_c / dH_h
            self._writeback_mflow_hot(mh_sol)
            return

        # ════════════════════════════════════════════════════════════════════
        # CASES: two enthalpy unknowns, cross-side
        # ════════════════════════════════════════════════════════════════════

        if unknowns == frozenset({'H_hi', 'H_ci'}):
            # free: T_hi; derive H_ci from energy balance
            def f(H_hi):
                try: H_hi = H_hi[0]
                except (TypeError, IndexError): pass
                try:
                    Q_h = mh * (H_hi - HO.H)
                    H_ci = CO.H - Q_h / mc
                    self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H=H_ci,P = CI.P,mflow= mc)
                    return _res_cth(Q_h, mh, mc)
                except Exception as e :
                    return _PENALTY * sign_factor

            Max_H_ci = self.Model.Prop(CI.fluid,StatePointName='_h_bound',P=CI.P, T = CP.PropsSI('TMAX', CI.fluid)).H
            Min_H_ci = self._min_H(CI, margin=0.0)
            H_lo = HO.H + mc*(CO.H - Max_H_ci)/mh+1e-3
            H_hi_bound = HO.H + mc*(CO.H - Min_H_ci)/mh-1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:                
                H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
            except Exception:                
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_ci={H_sol:.2f} J/kg")


            self.Hot_In = _h_from_H(mdl, HI.fluid, HI.StatePointName, H_sol, HI.P,mh)
            Q_h = mh * (self.Hot_In.H - HO.H)
            H_ci = CO.H - Q_h / mc
            self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H_ci, CI.P, mc)
            return
        # This cobination wont work. don't give reliable result.
        # if unknowns == frozenset({'H_hi', 'H_co'}):
        #     sign_factor = 1
        #     def f(H_hi):
        #         try: H_hi = H_hi[0] 
        #         except: pass
        #         try:
        #             Q_h = mh * (H_hi - HO.H)
        #             H_co = CI.H + Q_h / mc
        #             self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H=H_co,P = CO.P,mflow= mc)
        #             return _res_cth(Q_h, mh, mc)
        #         except Exception as e :
        #             print(e)
        #             return _PENALTY * sign_factor

        #     Max_H_co = self.Model.Prop(CO.fluid,StatePointName='_h_bound',P=CO.P, T = CP.PropsSI('TMAX', CO.fluid)).H
        #     Min_H_co = self.Model.Prop(CO.fluid,StatePointName='_l_bound',P=CO.P, T = CP.PropsSI('TMIN', CO.fluid)).H
        #     H_lo = HO.H + mc*(Min_H_co - CI.H)/mh+1e-3
        #     H_hi_bound = HO.H + mc*(Max_H_co - CI.H)/mh-1e-3
        #     Hi = f(H_hi_bound)
        #     Lo =  f(H_lo)
        #     print("AAAAAA",H_lo,H_hi_bound,Lo,Hi)
        #     if Lo == 1e10 and Hi > 0:
        #         sign_factor = -1  # try inverting the search direction if both ends invalid
        #     elif Lo == 1e10 and Hi < 0:
        #         sign_factor = 1
        #     elif Lo > 0 and Hi  == 1e10:
        #         sign_factor = -1
        #     elif Lo < 0 and Hi  == 1e10:
        #         sign_factor = 1
        #     try:                
        #         H_sol = brentq(f, H_lo, H_hi_bound, xtol=1e-4, rtol=1e-4)
        #     except:                
        #         H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
        #         print(f"Robust solver found H_ci={H_sol:.2f} K for {self.ID}")
        #     self.Hot_In = _h_from_H(mdl, HI.fluid, HI.StatePointName, H_sol, HI.P,mh)
        #     Q_h = mh * (self.Hot_In.H - HO.H)
        #     H_co = CI.H + Q_h / mc
        #     self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H_co, CO.P, mc)
        #     return

        # if unknowns == frozenset({'H_ho', 'H_ci'}):
        #     def f(T_ho):
        #         self.Hot_Out = _h_from_T(mdl, HO.fluid, HO.StatePointName, T_ho, HO.P)
        #         Q_h = mh * (HI.H - self.Hot_Out.H)
        #         H_ci = CO.H - Q_h / mc
        #         self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H_ci, CI.P)
        #         return _res_cth(Q_h, mh, mc)
        #     T_lo = CP.PropsSI('TMIN', HO.fluid) + 1
        #     T_sol = brentq(f, T_lo, HI.T - 1e-3, xtol=1e-4, rtol=1e-4)
        #     self.Hot_Out = _h_from_T(mdl, HO.fluid, HO.StatePointName, T_sol, HO.P,mh)
        #     Q_h = mh * (HI.H - self.Hot_Out.H)
        #     H_ci = CO.H - Q_h / mc
        #     self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H_ci, CI.P, mc)
        #     return

        if unknowns == frozenset({'H_ho', 'H_co'}):
            def f(H_ho):
                try:
                    self.Hot_Out = _h_from_H(mdl, HO.fluid, HO.StatePointName, H_ho, HO.P,mh)
                    Q_h = mh * (HI.H - self.Hot_Out.H)
                    H_co = CI.H + Q_h / mc
                    self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H_co, CO.P,mc)
                    return _res_cth(Q_h, mh, mc)
                except Exception as e:
                    return _PENALTY * sign_factor
            Max_H_co = self.Model.Prop(CO.fluid,StatePointName='_h_bound',P=CO.P, T = CP.PropsSI('TMAX', CO.fluid)).H
            Min_H_co = self._min_H(CO, margin=0.0)
            H_lo = HI.H - mc*(Max_H_co - CI.H )/mh + 1e-3
            H_hi_bound = HI.H - mc*(Min_H_co - CI.H )/mh - 1e-3
            sign_factor = _fun_sign_fac(f,H_lo,H_hi_bound)
            try:
                H_sol = brentq(f,H_lo,H_hi_bound , xtol=1e-4, rtol=1e-4)
            except Exception:
                H_sol = robust_solver(f, H_lo, H_hi_bound, tol=1e-4)
                if self.print_residual:
                    print(f"[{self.ID}] robust solver found H_ci={H_sol:.2f} J/kg")

            self.Hot_Out = _h_from_H(mdl, HO.fluid, HO.StatePointName, H_sol, HO.P,mh)
            Q_h = mh * (HI.H - self.Hot_Out.H)
            H_co = CI.H + Q_h / mc
            self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H_co, CO.P, mc)
            return

        raise ValueError(
            f"[{self.ID}] 2-unknown case not recognised: {unknowns}. "
            f"This combination may be unsolvable."
        )

    # ─────────────────────────────────────────────────────────────────────────
    #  3-UNKNOWN SOLVER ENGINE
    # ─────────────────────────────────────────────────────────────────────────

    def _anchor_rules(self):
        """
        Return list of (known_var, resolves_var, valid_hex_types) tuples
        describing all possible PPT diagonal anchors.
        """
        return [
            ('H_hi', 'H_co', ('double_pipe', 'Evaporator')),
            ('H_ho', 'H_ci', ('double_pipe', 'Condenser')),
            ('H_ci', 'H_ho', ('double_pipe', 'Condenser')),
            ('H_co', 'H_hi', ('double_pipe', 'Evaporator')),
        ]

    def _get_H(self, var_name):
        """Return enthalpy value for a named variable (None if unknown)."""
        mapping = {
            'H_hi': self.Hot_In.H  if self.Hot_In  else None,
            'H_ho': self.Hot_Out.H if self.Hot_Out else None,
            'H_ci': self.Cold_In.H if self.Cold_In else None,
            'H_co': self.Cold_Out.H if self.Cold_Out else None,
        }
        return mapping.get(var_name)

    def _apply_anchor(self, known_var, resolves_var):
        """
        Apply PPT diagonal anchor:
          known_var  → resolves_var
        Sets the resolved state point on self.

        Returns the anchor point index in the T-profile (0 = hot outlet end,
        N = hot inlet end) so the 2-unknown solver can exclude it.
        """
        mdl = self.Model
        HI  = self.Hot_In;  HO = self.Hot_Out
        CI  = self.Cold_In; CO = self.Cold_Out

        def _require_T(point, label):
            """Fail loudly (and usefully) instead of comparing None to a float."""
            if point is None or point.T is None:
                raise ValueError(
                    f"[{self.ID}] Cannot apply the {known_var} → {resolves_var} "
                    f"anchor: {label} temperature is unknown. The PPT diagonal "
                    f"needs a resolved temperature at the anchoring end."
                )
            return point.T

        if self.HEX_type in ('double_pipe', 'Condenser'):
            if known_var == 'H_ho' and resolves_var == 'H_ci':
                # Anchor T_ci = T_ho − PPT.  Only valid if the resulting cold
                # inlet stays BELOW the cold outlet, otherwise the HEX would
                # be driven to the unphysical state T_co < T_ci.
                T_ho = _require_T(HO, 'hot-outlet')
                T_co = _require_T(CO, 'cold-outlet')
                if T_ho - self.PPT >= T_co:
                    raise ValueError(
                        f"[{self.ID}] Cannot anchor H_ho → H_ci: it would set "
                        f"T_ci = {T_ho - self.PPT - 273.15:.2f} °C, which is not "
                        f"below T_co = {T_co - 273.15:.2f} °C."
                    )
                self.Cold_In = _h_from_T(mdl, CI.fluid, CI.StatePointName,
                                         T_ho - self.PPT, CI.P,
                                         self.Cold_Mass_flowrate)
                return 0   # anchor at hot-outlet end (index 0)

            if known_var == 'H_ci' and resolves_var == 'H_ho':
                T_ci = _require_T(CI, 'cold-inlet')
                self.Hot_Out = _h_from_T(mdl, HO.fluid, HO.StatePointName,
                                         T_ci + self.PPT, HO.P,
                                         self.Hot_Mass_flowrate)
                return 0   # anchor at hot-outlet end (index 0)

        if self.HEX_type in ('double_pipe', 'Evaporator'):
            if known_var == 'H_hi' and resolves_var == 'H_co':
                T_hi = _require_T(HI, 'hot-inlet')
                self.Cold_Out = _h_from_T(mdl, CO.fluid, CO.StatePointName,
                                          T_hi - self.PPT, CO.P,
                                          self.Cold_Mass_flowrate)
                return self.div_N   # anchor at hot-inlet end (index N)

            if known_var == 'H_co' and resolves_var == 'H_hi':
                # Anchor T_hi = T_co + PPT.  Only valid if the resulting hot
                # inlet stays ABOVE the hot outlet, otherwise the HEX would be
                # driven to the unphysical state T_hi < T_ho.
                # (v3.2.1 reused the *hot-outlet vs cold-outlet* test from the
                # branch above here, which tests the wrong pair entirely.)
                T_co = _require_T(CO, 'cold-outlet')
                T_ho = _require_T(HO, 'hot-outlet')
                if T_co + self.PPT <= T_ho:
                    raise ValueError(
                        f"[{self.ID}] Cannot anchor H_co → H_hi: it would set "
                        f"T_hi = {T_co + self.PPT - 273.15:.2f} °C, which is not "
                        f"above T_ho = {T_ho - 273.15:.2f} °C."
                    )
                self.Hot_In = _h_from_T(mdl, HI.fluid, HI.StatePointName,
                                        T_co + self.PPT, HI.P,
                                        self.Hot_Mass_flowrate)
                return self.div_N   # anchor at hot-inlet end (index N)

        raise ValueError(
            f"[{self.ID}] Anchor {known_var} → {resolves_var} is not defined "
            f"for HEX_type='{self.HEX_type}'."
        )

    def _snapshot(self):
        """Capture the four state points and both mass flows."""
        return (self.Hot_In, self.Hot_Out, self.Cold_In, self.Cold_Out,
                self.Hot_Mass_flowrate, self.Cold_Mass_flowrate)

    def _restore(self, snap):
        """Undo any mutation made while attempting a solve."""
        (self.Hot_In, self.Hot_Out, self.Cold_In, self.Cold_Out,
         self.Hot_Mass_flowrate, self.Cold_Mass_flowrate) = snap

    def _solve_3_unknown(self, unknowns, dP_h, dP_c):
        """
        3-unknown solver:
          1. Find a valid anchor for this HEX type
          2. Check that the anchor reduces the problem to a SOLVABLE
             2-unknown case *before* mutating anything
          3. Apply anchor → resolve 1 unknown
          4. Call _solve_2_unknown with exclude_idx set to the anchor point

        Step 2 is the important one.  Up to v3.2.1 the anchor was applied and
        written straight into ``Model.Point`` first; when the reduced problem
        turned out to be one of the unsolvable cases (2-unknown cases 8, 9 or
        15) the solver raised *after* the model had already been mutated,
        leaving the user with a silently corrupted state-point table.
        Nothing is committed here until the whole chain is known to work.
        """
        snap = self._snapshot()

        # Combinations that reduce to an unsolvable 2-unknown case.
        # See HEX_Logic_Table.xlsx, sheet "3-Unknown Cases": rows 1, 2, 5, 6
        # collapse to cases 9/8, rows 11-14 have no anchor at all, and rows
        # 17-20 collapse to case 15 (mh and mc hold a constant ratio, so
        # there are infinitely many solutions).
        candidates = []
        for known_var, resolves_var, valid_types in self._anchor_rules():
            if self.HEX_type not in valid_types:
                continue
            if known_var in unknowns:          # known_var must be known
                continue
            if resolves_var not in unknowns:   # resolves_var must be unknown
                continue
            if self._get_H(known_var) is None:
                continue
            candidates.append((known_var, resolves_var))

        if not candidates:
            raise ValueError(
                f"[{self.ID}] No valid anchor exists for unknowns="
                f"{set(unknowns)} with HEX_type='{self.HEX_type}'. "
                f"This combination is unsolvable: the PPT diagonal rule needs "
                f"one known enthalpy that resolves one of the unknowns."
            )

        # Keep only anchors that leave a solvable 2-unknown remainder.
        viable = [(k, r) for k, r in candidates
                  if (unknowns - {r}) in self._SOLVABLE_2_UNKNOWN]

        if not viable:
            reduced = {r: set(unknowns - {r}) for _, r in candidates}
            raise ValueError(
                f"[{self.ID}] Unsolvable 3-unknown combination "
                f"{set(unknowns)}: every available anchor reduces it to an "
                f"unsolvable 2-unknown case ({reduced}). No state points were "
                f"modified."
            )

        last_error = None
        for known_var, resolves_var in viable:
            try:
                anchor_idx = self._apply_anchor(known_var, resolves_var)
                remaining  = unknowns - {resolves_var}
                if self.print_residual:
                    print(f"[{self.ID}] anchored {known_var} → {resolves_var}; "
                          f"remaining: {set(remaining)}")
                self._solve_2_unknown(remaining, dP_h, dP_c,
                                      exclude_idx=anchor_idx)
                return
            except Exception as exc:
                # Roll back and try the next anchor (the logic table lists two
                # valid anchors for most solvable rows).
                self._restore(snap)
                last_error = exc

        raise ValueError(
            f"[{self.ID}] All {len(viable)} anchor(s) failed for unknowns="
            f"{set(unknowns)}. Last error: {last_error}"
        ) from last_error

    # ─────────────────────────────────────────────────────────────────────────
    #  CLASSIFY UNKNOWNS
    # ─────────────────────────────────────────────────────────────────────────

    def _identify_unknowns(self):
        """Return frozenset of unknown variable names."""
        unknowns = set()
        if self.Hot_In  is None or self.Hot_In.H  is None: unknowns.add('H_hi')
        if self.Hot_Out is None or self.Hot_Out.H is None: unknowns.add('H_ho')
        if self.Cold_In  is None or self.Cold_In.H  is None: unknowns.add('H_ci')
        if self.Cold_Out is None or self.Cold_Out.H is None: unknowns.add('H_co')
        if self.Hot_Mass_flowrate  is None: unknowns.add('mh')
        if self.Cold_Mass_flowrate is None: unknowns.add('mc')
        return frozenset(unknowns)

    # ─────────────────────────────────────────────────────────────────────────
    #  MAIN CALCULATION
    # ─────────────────────────────────────────────────────────────────────────

    def Cal(self):
        # Fetch state points
        self.Hot_In   = self.Model.Point[self.Hot_In_state]   if self.Hot_In_state   else None
        self.Hot_Out  = self.Model.Point[self.Hot_Out_state]  if self.Hot_Out_state  else None
        self.Cold_In  = self.Model.Point[self.Cold_In_state]  if self.Cold_In_state  else None
        self.Cold_Out = self.Model.Point[self.Cold_Out_state] if self.Cold_Out_state else None
        
        # Resolve mass flows from state points.
        # This MUST happen before any solver branch, including SimpleHEX:
        # up to v3.2.1 the SimpleHEX branch returned above this block, so
        # Hot/Cold_Mass_flowrate were still None from __init__ and the
        # solver silently took the wrong branch.
        if self.Hot_In is not None and self.Hot_Out is not None:
            self.Hot_Mass_flowrate = self._resolve_mass_flow(
                self.Hot_In, self.Hot_Out, 'Hot')
        if self.Cold_In is not None and self.Cold_Out is not None:
            self.Cold_Mass_flowrate = self._resolve_mass_flow(
                self.Cold_In, self.Cold_Out, 'Cold')

        # SimpleHEX branch
        if self.HEX_type == 'SimpleHEX':
            self._solve_simple_hex()
            self._compute_exergy_destruction()
            self.Solution_Status = True
            self._writeback_all()
            return

        self._check_PPT_ends()

        # Pressure drops
        dP_h = (self.Hot_In.P  - self.Hot_Out.P)  if (self.Hot_In  and self.Hot_Out)  else 0.0
        dP_c = (self.Cold_In.P - self.Cold_Out.P) if (self.Cold_In and self.Cold_Out) else 0.0

        # Identify unknowns
        unknowns = self._identify_unknowns()
        n = len(unknowns)
        if self.print_residual:
            print(f"[{self.ID}] Unknowns: {unknowns} (n={n})")
        # The residual functions handed to brentq assign to self.Hot_Out /
        # self.Cold_In / self.Cold_Out on EVERY probe, so a solver that raises
        # part-way leaves the component holding whatever the last (failed)
        # trial produced.  Snapshot up front and roll back on any failure, so
        # a caught exception never leaves a half-probed state behind for the
        # next Cal() or for the fallback path to build on.
        snap = self._snapshot()
        try:
            self._dispatch_solver(unknowns, n, dP_h, dP_c)

            # Post-solve validation and derived quantities
            self._check_energy_balance()
            self._check_PPT_ends()
            self._compute_outputs(dP_h, dP_c)
            self._compute_exergy_destruction()
        except Exception:
            self._restore(snap)
            raise

        self.Solution_Status = True
        self._writeback_all()

    def _dispatch_solver(self, unknowns, n, dP_h, dP_c):
        """Route to the solver that matches the number of unknowns."""
        # ── NTU-ε path ──────────────────────────────────────────────────
        # This MUST be tested before the unknown-count dispatch.  In v3.2.1
        # it was the last `elif` in a chain that had already consumed
        # n == 0, 1, 2 and 3, so `_solve_effectiveness()` was unreachable and
        # the `effectiveness` constructor argument silently did nothing.
        if (self.effectiveness is not None
                and unknowns <= frozenset({'H_ho', 'H_co'})):
            self._solve_effectiveness()

        elif n == 0:
            # All known — consistency check + fill profile
            self._check_energy_balance()
            Q = (self.Hot_In.H - self.Hot_Out.H) * self.Hot_Mass_flowrate
            Th, Tc, dT, _ = self._build_profile(
                Q, dP_h, dP_c,
                self.Hot_Out.H, self.Cold_In.H,
                self.Hot_Out.P, self.Cold_In.P,
                self.Hot_Mass_flowrate, self.Cold_Mass_flowrate,
                'cold_to_hot')

        elif n == 1:
            unk = next(iter(unknowns))
            self._solve_1_unknown(unk, dP_h, dP_c)

        elif n == 2:
            self._solve_2_unknown(unknowns, dP_h, dP_c)

        elif n == 3:
            self._solve_3_unknown(unknowns, dP_h, dP_c)

        else:
            raise ValueError(
                f"[{self.ID}] Underdetermined: {n} unknowns {set(unknowns)}. "
                f"Maximum supported is 3. Supply more inputs."
            )

    # ─────────────────────────────────────────────────────────────────────────
    #  1-UNKNOWN SOLVER  (direct / single brentq)
    # ─────────────────────────────────────────────────────────────────────────

    def _solve_1_unknown(self, unk, dP_h, dP_c):
        mdl = self.Model
        HI = self.Hot_In;  HO = self.Hot_Out
        CI = self.Cold_In; CO = self.Cold_Out
        mh = self.Hot_Mass_flowrate
        mc = self.Cold_Mass_flowrate

        if unk == 'H_ho':
            Q_c = mc * (CO.H - CI.H)
            H_ho = HI.H - Q_c / mh
            self.Hot_Out = _h_from_H(mdl, HO.fluid, HO.StatePointName, H_ho, HO.P, mh)

        elif unk == 'H_hi':
            Q_c = mc * (CO.H - CI.H)
            H_hi = HO.H + Q_c / mh
            self.Hot_In = _h_from_H(mdl, HI.fluid, HI.StatePointName, H_hi, HI.P, mh)

        elif unk == 'H_co':
            Q_h = mh * (HI.H - HO.H)
            H_co = CI.H + Q_h / mc
            self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, H_co, CO.P, mc)

        elif unk == 'H_ci':
            Q_h = mh * (HI.H - HO.H)
            H_ci = CO.H - Q_h / mc
            self.Cold_In = _h_from_H(mdl, CI.fluid, CI.StatePointName, H_ci, CI.P, mc)

        elif unk == 'mh':
            Q_c = mc * (CO.H - CI.H)
            dH_h = HI.H - HO.H
            if abs(dH_h) < 1e-9:
                raise ValueError(f"[{self.ID}] Cannot solve mh: ΔH_hot = 0.")
            mh_sol = Q_c / dH_h
            self._writeback_mflow_hot(mh_sol)

        elif unk == 'mc':
            Q_h = mh * (HI.H - HO.H)
            dH_c = CO.H - CI.H
            if abs(dH_c) < 1e-9:
                raise ValueError(f"[{self.ID}] Cannot solve mc: ΔH_cold = 0.")
            mc_sol = Q_h / dH_c
            self._writeback_mflow_cold(mc_sol)

    # ─────────────────────────────────────────────────────────────────────────
    #  POST-SOLVE OUTPUTS
    # ─────────────────────────────────────────────────────────────────────────

    def _compute_outputs(self, dP_h, dP_c):
        HI = self.Hot_In;  HO = self.Hot_Out
        CI = self.Cold_In; CO = self.Cold_Out
        mh = self.Hot_Mass_flowrate
        mc = self.Cold_Mass_flowrate

        # Q
        # NOTE: these were truthiness tests up to v3.2.1 (`if mh and HI.H`).
        # CoolProp's reference state puts h ≈ 0 near the triple point, so a
        # perfectly valid H_hi = 0.0 read as "missing" and the duty was
        # silently skipped.  Likewise a legitimate zero-duty result (Q = 0)
        # suppressed the UA calculation below.
        if (mh is not None and HI is not None and HO is not None
                and HI.H is not None and HO.H is not None):
            self.Q = (HI.H - HO.H) * mh
        elif (mc is not None and CI is not None and CO is not None
                and CI.H is not None and CO.H is not None):
            self.Q = (CO.H - CI.H) * mc

        # Hot-to-cold mass flow ratio
        if (HI is not None and HO is not None
                and CI is not None and CO is not None
                and None not in (HI.H, HO.H, CI.H, CO.H)):
            dH_h = HI.H - HO.H
            dH_c = CO.H - CI.H
            if abs(dH_h) > 1e-9:
                self.Hot_to_Cold = dH_c / dH_h

        # UA via LMTD
        if (self.Q is not None and None not in (HI, HO, CI, CO)
                and None not in (HI.T, HO.T, CI.T, CO.T)):
            dT1 = HI.T  - CO.T
            dT2 = HO.T  - CI.T
            if dT1 > 0 and dT2 > 0:
                LMTD = ((dT1 - dT2) / np.log(dT1 / dT2)
                        if abs(dT1 - dT2) > 1e-6 else dT1)
                # LMTD was computed and discarded up to v3.2.1, even though
                # __str__ reports UA and users routinely want the mean
                # temperature difference that produced it.
                self.LMTD = LMTD
                if abs(LMTD) > 1e-9:
                    self.UA = self.Q / LMTD

        # Pinch graph
        if (self.PPT_graph and self.Q is not None
                and mh is not None and mc is not None):
            Th, Tc, dT, _ = self._build_profile(
                self.Q, dP_h, dP_c,
                HO.H, CI.H, HO.P, CI.P,
                mh, mc, 'cold_to_hot')
            
            self._draw_pinch(Th, Tc, dT)

    def _compute_exergy_destruction(self):
        """
        Exergy destruction, Ex_D = ΣEx_in − ΣEx_out − Ex_Q.

        Three cases, because a SimpleHEX may have only one connected side:

        * **Both sides connected.** Adiabatic to the surroundings, so
          Ex_Q = 0 and Ex_D = (Ex_hi + Ex_ci) − (Ex_ho + Ex_co).

        * **Hot side only** (``HeatAdded=False`` — a condenser or cooler
          rejecting heat to the environment). The sink is the dead state, so
          the rejected heat carries no exergy (Ex_Q = 0) and
          Ex_D = Ex_hi − Ex_ho.

        * **Cold side only** (``HeatAdded=True`` — a boiler or heater driven
          by an unmodelled source). Ex_D = Ex_Q,source − (Ex_co − Ex_ci), and
          Ex_Q,source depends on the source temperature, which this component
          does not know. Reporting −(Ex_co − Ex_ci) would be a *negative*
          "destruction" that poisons the total in ModelSummary, so the result
          is left as "Not Calculated". Model the source explicitly (as the hot
          side of a real HEX, or via a Source component) to get a number here.
        """
        hot  = self.Hot_In  is not None and self.Hot_Out  is not None
        cold = self.Cold_In is not None and self.Cold_Out is not None

        if cold and not hot:
            self.Ex_D = "Not Calculated"
            return

        legs = []
        if hot:
            legs.append((self.Hot_In, self.Hot_Out))
        if cold:
            legs.append((self.Cold_In, self.Cold_Out))

        if not legs:
            self.Ex_D = "Not Calculated"
            return

        total = 0.0
        for pt_in, pt_out in legs:
            if pt_in.Ex is None or pt_out.Ex is None:
                self.Ex_D = "Not Calculated"
                return
            total += pt_in.Ex - pt_out.Ex

        self.Ex_D = total

    def _writeback_all(self):
        if self.Hot_In_state:   self.Model.Point[self.Hot_In_state]   = self.Hot_In
        if self.Hot_Out_state:  self.Model.Point[self.Hot_Out_state]  = self.Hot_Out
        if self.Cold_In_state:  self.Model.Point[self.Cold_In_state]  = self.Cold_In
        if self.Cold_Out_state: self.Model.Point[self.Cold_Out_state] = self.Cold_Out

    # ─────────────────────────────────────────────────────────────────────────
    #  EFFECTIVENESS SOLVER  (NTU-ε)
    # ─────────────────────────────────────────────────────────────────────────

    def _solve_effectiveness(self):
        mdl = self.Model
        HI = self.Hot_In;  HO = self.Hot_Out
        CI = self.Cold_In; CO = self.Cold_Out
        mh = self.Hot_Mass_flowrate
        mc = self.Cold_Mass_flowrate

        # Effectiveness is defined as eps = Q / Q_max, where Q_max is the duty
        # of an infinitely long counter-flow exchanger: the limiting stream
        # leaves at the OTHER stream's inlet temperature.  In enthalpy terms
        #
        #     Q_max = min[ mh (h_hi - h_h(T_ci)),  mc (h_c(T_hi) - h_ci) ]
        #
        # The familiar Cmin (T_hi - T_ci) form is this expression with each
        # enthalpy difference replaced by cp dT, i.e. it assumes cp is
        # constant across the exchanger.  That assumption is convenient by
        # hand but unnecessary here, and it is poor for gases over a wide
        # temperature span: cp for air rises about 19 % between 450 K and
        # 1140 K, which is a routine regenerator duty.
        #
        # The enthalpy form is used unconditionally.  It is exact for real
        # fluids, handles two-phase and incompressible streams without a
        # special case, and reduces to Cmin (T_hi - T_ci) whenever cp really
        # is constant.
        #
        # Properties are taken through Model.Prop rather than CoolProp
        # directly so that the incompressible correlations (Therminol-66)
        # are handled by the same path as everything else.
        T_hi, T_ci = HI.T, CI.T
        h_hi_at_Tci = mdl.Prop(HI.fluid, StatePointName='_eps_h',
                               T=T_ci, P=HO.P).H
        h_ci_at_Thi = mdl.Prop(CI.fluid, StatePointName='_eps_c',
                               T=T_hi, P=CO.P).H
        Q_max = min((HI.H - h_hi_at_Tci) * mh,
                    (h_ci_at_Thi - CI.H) * mc)
        if Q_max < 0:
            raise ValueError(
                f"[{self.ID}] Effectiveness solver: the hot inlet "
                f"({T_hi - 273.15:.2f} °C) is not hotter than the cold inlet "
                f"({T_ci - 273.15:.2f} °C), so no heat can be transferred."
            )
        self.Q = self.effectiveness * Q_max

        hh = HI.H - self.Q / mh
        hc = CI.H + self.Q / mc
        self.Hot_Out  = _h_from_H(mdl, HO.fluid, HO.StatePointName, hh, HO.P, mh)
        self.Cold_Out = _h_from_H(mdl, CO.fluid, CO.StatePointName, hc, CO.P, mc)

    # ─────────────────────────────────────────────────────────────────────────
    #  SIMPLEHEX SOLVER  (carried over from v2.0)
    # ─────────────────────────────────────────────────────────────────────────

    def _solve_simple_hex(self):
        """
        SimpleHEX: a pure energy-balance heat exchanger with no pinch
        analysis and no temperature profile.  Either side may be omitted
        (``Hot_*_state=None`` or ``Cold_*_state=None``) to model a boiler,
        condenser, heater or cooler against an unmodelled utility.

        Each side is solved independently against the shared duty ``Q``.
        For a given side, exactly one of {inlet H, outlet H, mass flow, Q}
        may be unknown; everything else must be supplied.

        This replaces the single if/elif chain used up to v3.2.1, which had
        three defects:
          * only ONE branch could ever run, so a two-sided HEX with a known
            cold side and an unknown hot outlet computed Q and returned
            without ever touching the hot side;
          * a fully-determined side hit no branch at all and fell through to
            "insufficient inputs", so calling Cal() twice on a solved HEX
            raised — which made every iteration of Solve() report phantom
            failures;
          * the pre-computed self.Q from a previous call was never refreshed.
        """
        sides = []
        if self.Cold_In is not None and self.Cold_Out is not None:
            sides.append('cold')
        if self.Hot_In is not None and self.Hot_Out is not None:
            sides.append('hot')

        if not sides:
            raise ValueError(
                f"SimpleHEX {self.ID}: neither side is fully connected. "
                f"Provide both Hot_In_state/Hot_Out_state or both "
                f"Cold_In_state/Cold_Out_state (or both pairs)."
            )

        # ── Pass 1: derive Q from whichever side is fully determined ──────
        Q_from = {}
        for side in sides:
            q = self._simple_side_duty(side)
            if q is not None:
                Q_from[side] = q

        if Q_from:
            values = list(Q_from.values())
            if len(values) == 2:
                ref = max(abs(values[0]), abs(values[1]), 1.0)
                if abs(values[0] - values[1]) / ref > 0.01:
                    raise ValueError(
                        f"SimpleHEX {self.ID}: energy imbalance — "
                        f"Q_{list(Q_from)[0]}={values[0]:.2f} W vs "
                        f"Q_{list(Q_from)[1]}={values[1]:.2f} W "
                        f"({abs(values[0] - values[1]) / ref * 100:.2f} % apart)."
                    )
            Q_new = sum(values) / len(values)

            if self.Q is not None and len(values) == len(sides):
                # Over-specified: the user supplied Q *and* both endpoints.
                ref = max(abs(self.Q), abs(Q_new), 1.0)
                if abs(self.Q - Q_new) / ref > 0.01:
                    raise ValueError(
                        f"SimpleHEX {self.ID}: supplied Q={self.Q:.2f} W "
                        f"is inconsistent with the state points "
                        f"(Q={Q_new:.2f} W from the energy balance)."
                    )
            self.Q = Q_new

        if self.Q is None:
            raise ValueError(
                f"SimpleHEX {self.ID}: cannot determine the heat duty. "
                f"Supply Q explicitly, or fully define at least one side "
                f"(both enthalpies and the mass flow rate)."
            )

        # ── Pass 2: close every side that is still incomplete ─────────────
        for side in sides:
            self._simple_close_side(side)

    def _simple_side_duty(self, side):
        """Return Q for `side` if that side is fully determined, else None."""
        if side == 'cold':
            pt_in, pt_out, m = self.Cold_In, self.Cold_Out, self.Cold_Mass_flowrate
            sign = 1.0            # cold side gains energy
        else:
            pt_in, pt_out, m = self.Hot_In, self.Hot_Out, self.Hot_Mass_flowrate
            sign = -1.0           # hot side loses energy

        if None in (pt_in.H, pt_out.H, m):
            return None
        return sign * (pt_out.H - pt_in.H) * m

    def _simple_close_side(self, side):
        """Solve whichever single unknown remains on `side`, given self.Q."""
        mdl = self.Model
        if side == 'cold':
            pt_in, pt_out = self.Cold_In, self.Cold_Out
            m    = self.Cold_Mass_flowrate
            sign = 1.0
        else:
            pt_in, pt_out = self.Hot_In, self.Hot_Out
            m    = self.Hot_Mass_flowrate
            sign = -1.0

        unknowns = [n for n, v in (('H_in', pt_in.H), ('H_out', pt_out.H),
                                   ('m', m)) if v is None]

        if not unknowns:
            return                                  # already determined

        if len(unknowns) > 1:
            raise ValueError(
                f"SimpleHEX {self.ID}: the {side} side has "
                f"{len(unknowns)} unknowns {unknowns}. At most one of "
                f"(inlet H, outlet H, mass flow) may be unknown once Q is known."
            )

        unk = unknowns[0]

        if unk == 'H_out':
            H = pt_in.H + sign * self.Q / m
            new = _h_from_H(mdl, pt_out.fluid, pt_out.StatePointName,
                            H, pt_out.P, m)
            setattr(self, 'Cold_Out' if side == 'cold' else 'Hot_Out', new)

        elif unk == 'H_in':
            H = pt_out.H - sign * self.Q / m
            new = _h_from_H(mdl, pt_in.fluid, pt_in.StatePointName,
                            H, pt_in.P, m)
            setattr(self, 'Cold_In' if side == 'cold' else 'Hot_In', new)

        else:   # unk == 'm'
            dH = sign * (pt_out.H - pt_in.H)
            if abs(dH) < 1e-9:
                raise ValueError(
                    f"SimpleHEX {self.ID}: {side}-side ΔH is zero, so the "
                    f"mass flow rate cannot be derived from Q."
                )
            m_new = self.Q / dH
            if m_new <= 0:
                raise ValueError(
                    f"SimpleHEX {self.ID}: derived {side}-side mass flow is "
                    f"non-positive ({m_new:.6g} kg/s). Check the sign of Q "
                    f"and the HeatAdded flag."
                )
            if side == 'cold':
                self._writeback_mflow_cold(m_new)
            else:
                self._writeback_mflow_hot(m_new)

    # ─────────────────────────────────────────────────────────────────────────
    #  PINCH PLOT
    # ─────────────────────────────────────────────────────────────────────────

    def _draw_pinch(self, Th, Tc, dT):
        plt.figure()
        plt.plot(range(len(Th)), Th - 273.15, 'r-', label='Hot')
        plt.plot(range(len(Tc)), Tc - 273.15, 'b-', label='Cold')
        plt.axhline(y=np.min(Th - 273.15), color='gray', linestyle='--', alpha=0.4)
        y1,y2 = (min(min(Tc),min(Th))-273.15-10,max(max(Tc),max(Th))-273.15+10)
        plt.arrow(np.argmin(dT),(Tc[np.argmin(dT)]+Th[np.argmin(dT)])/2-273.15,-10,-(Tc[np.argmin(dT)]-273.15-y1)*0.6,head_width=2, head_length=5)
        plt.text(np.argmin(dT)-40, y1+1.5, f"PPT = {min(dT):.2f}\n", fontsize=12, color='blue')
        plt.ylim((y1,y2))
        plt.legend()
        plt.title(f"{self.ID}  —  min ΔT = {np.min(dT):.2f} K  |  PPT = {self.PPT:.2f} K")
        plt.xlabel('Segment')
        plt.ylabel('Temperature [°C]')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    # ─────────────────────────────────────────────────────────────────────────
    #  __str__
    # ─────────────────────────────────────────────────────────────────────────

    def __str__(self):
        if self.HEX_type == 'SimpleHEX':
            return (f"{self.ID} (SimpleHEX):\n"
                    f"  Q     : {self.Q} W\n"
                    f"  Solved: {self.Solution_Status}")
        try:
            return (
                f"{self.ID} ({self.HEX_type}):\n"
                f"  Hot  : {self.Hot_In.T-273.15:.1f} → {self.Hot_Out.T-273.15:.1f} °C"
                f"  mh={self.Hot_Mass_flowrate:.3f} kg/s\n"
                f"  Cold : {self.Cold_In.T-273.15:.1f} → {self.Cold_Out.T-273.15:.1f} °C"
                f"  mc={self.Cold_Mass_flowrate:.3f} kg/s\n"
                f"  Q    : {self.Q:.2f} W\n"
                f"  UA   : {self.UA} W/K\n"
                f"  LMTD : {self.LMTD} K\n"
                f"  Ex_D : {self.Ex_D}\n"
                f"  Solved: {self.Solution_Status}"
            )
        except Exception:
            return f"{self.ID}: solved={self.Solution_Status}"


# ─────────────────────────────────────────────────────────────────────────────
#  TES  (Thermal Energy Storage) — unchanged from v2.0
# ─────────────────────────────────────────────────────────────────────────────
class TES(Component):

    def __init__(self, Model, ID, PPT, Charge, T_melt,
                 Hot_In_state, Hot_Out_state,
                 Cold_In_state, Cold_Out_state,
                 Charging_time, Discharging_time,
                 per_loss, Capacity=None, Calculate=False):
        self.T_melt            = T_melt
        self.per_loss          = per_loss
        self.Hot_In_state      = Hot_In_state
        self.Hot_Out_state     = Hot_Out_state
        self.Cold_In_state     = Cold_In_state
        self.Cold_Out_state    = Cold_Out_state
        self.Charging_time     = Charging_time
        self.Discharging_time  = Discharging_time
        self.Charge            = Charge
        self.Capacity          = Capacity
        self.PPT               = PPT
        self.Charging_Power    = None
        self.Discharging_Power = None
        self.Hot_Mass_flowrate  = None
        self.Cold_Mass_flowrate = None
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.Hot_In   = self.Model.Point[self.Hot_In_state]   if self.Hot_In_state   else None
        self.Hot_Out  = self.Model.Point[self.Hot_Out_state]  if self.Hot_Out_state  else None
        self.Cold_In  = self.Model.Point[self.Cold_In_state]  if self.Cold_In_state  else None
        self.Cold_Out = self.Model.Point[self.Cold_Out_state] if self.Cold_Out_state else None

        self.Hot_Mass_flowrate  = None
        self.Cold_Mass_flowrate = None

        # Reuse HeatExchanger's resolver rather than a third hand-rolled copy.
        # The old inline version compared floats with `mi == mo`, so two flows
        # that differed by one ULP raised a spurious mismatch error.
        if self.Hot_In is not None and self.Hot_Out is not None:
            self.Hot_Mass_flowrate = HeatExchanger._resolve_mass_flow(
                self, self.Hot_In, self.Hot_Out, 'Hot')
        if self.Cold_In is not None and self.Cold_Out is not None:
            self.Cold_Mass_flowrate = HeatExchanger._resolve_mass_flow(
                self, self.Cold_In, self.Cold_Out, 'Cold')

        if self.Charge == 'Discharging':
            if self.Capacity is not None:
                self.CapacityD = self.Capacity * (1 - self.per_loss)
            else:
                self.CapacityD = None

            CI = self.Cold_In;  CO = self.Cold_Out
            mc = self.Cold_Mass_flowrate
            dt = self.Discharging_time * 3600

            if CI.H is not None and CO.H is None and self.CapacityD is None:
                if self.PPT < (self.T_melt - CI.T):
                    self.Cold_Out = self.Model.Prop(
                        CO.fluid, StatePointName=CO.StatePointName,
                        P=CO.P, T=self.T_melt - self.PPT, Mass_flowrate=mc)
                    self.CapacityD = (self.Cold_Out.H - CI.H) * mc * dt
                else:
                    raise ValueError(f"[{self.ID}] Discharging PPT violation.")

            elif CI.H is not None and CO.H is None and self.CapacityD is not None:
                H = CI.H + self.CapacityD / (mc * dt)
                self.Cold_Out = self.Model.Prop(
                    CO.fluid, StatePointName=CO.StatePointName, P=CO.P, H=H, Mass_flowrate=mc)

            elif CI.H is not None and CO.H is not None and self.CapacityD is None:
                self.CapacityD = (CO.H - CI.H) * mc * dt

            elif CI.H is not None and CO.H is not None and self.CapacityD is not None:
                Q_check = (CO.H - CI.H) * mc * dt
                tol = 0.01 * max(abs(self.CapacityD), 1.0)
                if abs(Q_check - self.CapacityD) > tol:
                    raise ValueError(f"[{self.ID}] TES Discharging over-specified and inconsistent.")

            elif CI.H is None and CO.H is None and self.CapacityD is not None:
                self.Cold_Out = self.Model.Prop(
                    CO.fluid, StatePointName=CO.StatePointName,
                    P=CO.P, T=self.T_melt - self.PPT, Mass_flowrate=mc)
                H = self.Cold_Out.H - self.CapacityD / (mc * dt)
                self.Cold_In = self.Model.Prop(
                    CI.fluid, StatePointName=CI.StatePointName, P=CI.P, H=H, Mass_flowrate=mc)

            elif CI.H is None and CO.H is not None and self.CapacityD is not None:
                H = CO.H - self.CapacityD / (mc * dt)
                self.Cold_In = self.Model.Prop(
                    CI.fluid, StatePointName=CI.StatePointName, P=CI.P, H=H, Mass_flowrate=mc)
            else:
                raise ValueError(f"[{self.ID}] TES Discharging: unresolvable combination.")

            self.Discharging_Power = self.CapacityD / dt

        elif self.Charge == 'Charging':
            HI = self.Hot_In;  HO = self.Hot_Out
            mh = self.Hot_Mass_flowrate
            dt = self.Charging_time * 3600

            if HI.H is not None and HO.H is None and self.Capacity is None:
                if self.PPT < (HI.T - self.T_melt):
                    self.Hot_Out = self.Model.Prop(
                        HO.fluid, StatePointName=HO.StatePointName,
                        P=HO.P, T=self.T_melt + self.PPT, Mass_flowrate=mh)
                    self.Capacity = (HI.H - self.Hot_Out.H) * mh * dt
                else:
                    raise ValueError(f"[{self.ID}] Charging PPT violation.")

            elif HI.H is not None and HO.H is None and self.Capacity is not None:
                H = HI.H - self.Capacity / (mh * dt)
                self.Hot_Out = self.Model.Prop(
                    HO.fluid, StatePointName=HO.StatePointName, P=HO.P, H=H, Mass_flowrate=mh)

            elif HI.H is not None and HO.H is not None and self.Capacity is None:
                self.Capacity = (HI.H - HO.H) * mh * dt

            elif HI.H is not None and HO.H is not None and self.Capacity is not None:
                Q_check = (HI.H - HO.H) * mh * dt
                tol = 0.01 * max(abs(self.Capacity), 1.0)
                if abs(Q_check - self.Capacity) > tol:
                    raise ValueError(f"[{self.ID}] TES Charging over-specified and inconsistent.")

            elif HI.H is None and HO.H is None and self.Capacity is not None:
                self.Hot_Out = self.Model.Prop(
                    HO.fluid, StatePointName=HO.StatePointName,
                    P=HO.P, T=self.T_melt + self.PPT, Mass_flowrate=mh)
                H = self.Hot_Out.H + self.Capacity / (mh * dt)
                self.Hot_In = self.Model.Prop(
                    HI.fluid, StatePointName=HI.StatePointName, P=HI.P, H=H, Mass_flowrate=mh)

            elif HI.H is None and HO.H is not None and self.Capacity is not None:
                H = HO.H + self.Capacity / (mh * dt)
                self.Hot_In = self.Model.Prop(
                    HI.fluid, StatePointName=HI.StatePointName, P=HI.P, H=H, Mass_flowrate=mh)
            else:
                raise ValueError(f"[{self.ID}] TES Charging: unresolvable combination.")

            self.Charging_Power = self.Capacity / dt
        else:
            raise ValueError(f"Invalid Charge='{self.Charge}' in {self.ID}.")

        try:
            self.Ex_D = ((self.Hot_In.Ex + self.Cold_In.Ex)
                         - (self.Hot_Out.Ex + self.Cold_Out.Ex))
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True
        if self.Hot_In_state:   self.Model.Point[self.Hot_In_state]   = self.Hot_In
        if self.Hot_Out_state:  self.Model.Point[self.Hot_Out_state]  = self.Hot_Out
        if self.Cold_In_state:  self.Model.Point[self.Cold_In_state]  = self.Cold_In
        if self.Cold_Out_state: self.Model.Point[self.Cold_Out_state] = self.Cold_Out

    def __str__(self):
        return (
            f"{self.ID} (TES - {self.Charge}):\n"
            f"  Charging Power   : {self.Charging_Power} W\n"
            f"  Discharging Power: {self.Discharging_Power} W\n"
            f"  Capacity         : {self.Capacity} J\n"
            f"  Solved           : {self.Solution_Status}"
        )
