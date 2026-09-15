"""
ThermoSim validation: single-effect NH3-H2O absorption chiller.

Reproduces NH3H2O/demo_cycle.py's worked example -- same P_hi/P_lo, T_gen/T_abs
and rectified-vapor composition y_r -- but built entirely from ThermoSim's
Model/Prop/component API instead of hand-written balance equations, to
validate the mixture integration end to end:

  * Prop / NH3H2OFluid          -- the mixture state layer (state.py)
  * Pump, Expansion_valve       -- reused unchanged for the solution pump
                                    and both throttles
  * Mixer                       -- the absorber's adiabatic mixing step
  * Separator                   -- the generator's vapor/liquid split
                                    (mixture-aware fix, simple_components.py)
  * HeatExchanger(SimpleHEX)    -- generator pre-heat, condenser, evaporator
                                    and the absorber's cooling step ("Absorber
                                    = Mixer + HeatExchanger acting like a
                                    condenser", "Evaporator = HeatExchanger")

Rectifier
---------
Not modelled as a component. Like demo_cycle.py itself, the vapor leaving
the generator is assumed purified to a fixed composition y_r = 0.995 by an
unmodelled rectifying column, and the resulting vapor/weak-liquid mass split
follows from an overall NH3 species balance across the combined
generator+rectifier control volume (`f`, below) -- not from the raw
generator flash's own quality. That combined control volume has one inlet
splitting into two differently-flowing outlets, which does not fit the
HeatExchanger's single-stream-per-side port model, so its total duty
(`Q_genset`) is computed directly from the already-solved state points,
exactly as demo_cycle.py computes it -- not through a HeatExchanger call.
The physical generator HEX + Separator pair is still used to obtain the
weak-liquid composition x_w (= the generator flash's true liquid-phase
composition), which is the one quantity the rectifier idealization assumes
is unaffected by reflux.

Reference values (from running NH3H2O/demo_cycle.py directly):
  Q_gen (gen+rect) =   723.5428 kW      Q_cond =   327.9912 kW
  Q_evap           =   322.4703 kW      Q_abs  =   719.3540 kW
  W_pump           =     1.3322 kW      COP    =     0.4449
"""
import warnings
warnings.filterwarnings('ignore')

from ThermoSim import (
    ThermodynamicModel, NH3H2OFluid,
    Pump, Expansion_valve, Mixer, Separator, HeatExchanger,
)
import ThermoSim.nh3h2o as aw
from scipy.optimize import brentq

# ── design point (matches NH3H2O/demo_cycle.py) ────────────────────────────
P_hi, P_lo = 15.0e5, 4.0e5      # Pa
T_gen, T_abs = 400.0, 310.0     # K
y_r = 0.995                     # NH3 mass fraction leaving the rectifier
m_s = 1.0                       # kg/s of strong solution

# Absorber design composition: the z whose bubble point at P_lo is T_abs.
z_s = brentq(lambda z: aw.bubble_T(P_lo / 1e5, z) - T_abs, 0.01, 0.99, xtol=1e-12)

M = ThermodynamicModel()

# ── 1: absorber outlet -> 2: pump outlet ────────────────────────────────────
st1 = M.add_point(NH3H2OFluid(z_s), 'st1', P=P_lo, Q=0.0, Mass_flowrate=m_s)
st2 = M.add_point(NH3H2OFluid(z_s), 'st2', P=P_hi)
Pump(M, 'SolutionPump', In_state='st1', Out_state='st2',
     n_isen=1.0, Calculate=True)
st2 = M.Point['st2']   # Pump.Cal() replaces the Prop object; re-fetch it

# ── generator: heat the strong solution to T_gen, then split ───────────────
gen = M.add_point(NH3H2OFluid(z_s), 'gen', P=P_hi, T=T_gen, Mass_flowrate=m_s)
HeatExchanger(M, 'Generator', PPT=None, HEX_type='SimpleHEX', HeatAdded=True,
              Hot_In_state=None, Hot_Out_state=None,
              Cold_In_state='st2', Cold_Out_state='gen', Calculate=True)

gen_vap_raw = M.add_point(NH3H2OFluid(0.5), 'gen_vap_raw')   # placeholder; overwritten
gen_liq     = M.add_point(NH3H2OFluid(0.5), 'gen_liq')       # placeholder; overwritten
GenSep = Separator(M, 'GenSep', In_state='gen',
                    Out_vap_state='gen_vap_raw', Out_liq_state='gen_liq',
                    Calculate=True)
x_w = GenSep.Out_liq.z

# ── idealized rectifier: species balance fixes the vapor/liquid split ──────
f = (z_s - x_w) / (y_r - x_w)
m_v, m_w = f * m_s, (1.0 - f) * m_s
M.Point['gen_liq'].Mass_flowrate = m_w   # override the raw-flash split

rect_vap = M.add_point(NH3H2OFluid(y_r), 'rect_vap', P=P_hi, Q=1.0,
                        Mass_flowrate=m_v)

# ── condenser ────────────────────────────────────────────────────────────
cond_out = M.add_point(NH3H2OFluid(y_r), 'cond_out', P=P_hi, Q=0.0,
                        Mass_flowrate=m_v)
HeatExchanger(M, 'Condenser', PPT=None, HEX_type='SimpleHEX', HeatAdded=False,
              Hot_In_state='rect_vap', Hot_Out_state='cond_out',
              Cold_In_state=None, Cold_Out_state=None, Calculate=True)

# ── throttle 4->5, evaporator, throttle 7->8 ────────────────────────────────
ev_in = M.add_point(NH3H2OFluid(y_r), 'ev_in', P=P_lo, Mass_flowrate=m_v)
Expansion_valve(M, 'RefrigValve', In_state='cond_out', Out_state='ev_in',
                Calculate=True)

ev_out = M.add_point(NH3H2OFluid(y_r), 'ev_out', P=P_lo, Q=1.0,
                      Mass_flowrate=m_v)
HeatExchanger(M, 'Evaporator', PPT=None, HEX_type='SimpleHEX', HeatAdded=True,
              Hot_In_state=None, Hot_Out_state=None,
              Cold_In_state='ev_in', Cold_Out_state='ev_out', Calculate=True)

weak_lo = M.add_point(NH3H2OFluid(x_w), 'weak_lo', P=P_lo, Mass_flowrate=m_w)
Expansion_valve(M, 'SolutionValve', In_state='gen_liq', Out_state='weak_lo',
                Calculate=True)

# ── absorber: Mixer (adiabatic mix) + HeatExchanger acting as a condenser ──
mixed = M.add_point(NH3H2OFluid(z_s), 'mixed')
Mixer(M, 'Absorber_Mix', In_states=['ev_out', 'weak_lo'], Out_state='mixed',
      Calculate=True)
HeatExchanger(M, 'Absorber_Cool', PPT=None, HEX_type='SimpleHEX', HeatAdded=False,
              Hot_In_state='mixed', Hot_Out_state='st1',
              Cold_In_state=None, Cold_Out_state=None, Calculate=True)

# ── results ─────────────────────────────────────────────────────────────
Q_gen_raw = M.Component['Generator'].Q          # generator-only duty (real split)
Q_genset  = m_v * rect_vap.H + m_w * M.Point['gen_liq'].H - m_s * st2.H
Q_cond    = M.Component['Condenser'].Q
Q_evap    = M.Component['Evaporator'].Q
Q_abs     = M.Component['Absorber_Cool'].Q
W_pump    = m_s * (st2.H - st1.H)
COP       = Q_evap / (Q_genset + W_pump)

REF = {
    'Q_gen (gen+rect) kW': 723.5428, 'Q_cond kW': 327.9912,
    'Q_evap kW': 322.4703, 'Q_abs kW': 719.3540,
    'W_pump kW': 1.3322, 'COP': 0.4449,
}
GOT = {
    'Q_gen (gen+rect) kW': Q_genset / 1e3, 'Q_cond kW': Q_cond / 1e3,
    'Q_evap kW': Q_evap / 1e3, 'Q_abs kW': Q_abs / 1e3,
    'W_pump kW': W_pump / 1e3, 'COP': COP,
}

print('SINGLE-EFFECT NH3-H2O ABSORPTION CHILLER (built via ThermoSim components)')
print(f'  P_hi = {P_hi/1e5:g} bar, P_lo = {P_lo/1e5:g} bar, '
      f'T_gen = {T_gen:.0f} K, T_abs = {T_abs:.0f} K')
print(f'  strong z = {z_s:.5f}   weak x = {x_w:.5f}   refrigerant y = {y_r:.4f}')
print(f'  circulation ratio m_s/m_v = {1/f:.3f}')
print(f'  generator-only duty (real split, pre-rectifier) = {Q_gen_raw/1e3:.4f} kW')
print()
print(f'  {"quantity":<22s}{"ThermoSim":>12s}{"demo_cycle.py":>16s}{"diff %":>10s}')
for k, ref in REF.items():
    got = GOT[k]
    diff = 100 * abs(got - ref) / max(abs(ref), 1e-9)
    print(f'  {k:<22s}{got:12.4f}{ref:16.4f}{diff:10.4f}')

print()
tol_pct = 0.05
mismatches = [k for k, ref in REF.items()
              if 100 * abs(GOT[k] - ref) / max(abs(ref), 1e-9) > tol_pct]
if mismatches:
    raise AssertionError(f'Deviates from demo_cycle.py by more than '
                         f'{tol_pct}% for: {mismatches}')
print(f'All quantities match demo_cycle.py to within {tol_pct}%.')
