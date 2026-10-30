"""
Figure — organic Rankine cycle recovering heat from geothermal brine.

(a) Evaporator temperature profile at the design point.  The working-fluid
    mass flow rate is not supplied; it is returned by the exchanger from
    the specified minimum approach.
(b) Thermal efficiency and net power against evaporating pressure.

Regenerates fig_orc_example.png.

Run with:  MPLBACKEND=Agg python fig_orc_example.py
"""
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from ThermoSim import ThermodynamicModel, Turbine, Pump, HeatExchanger
from ThermoSim.plotting import CyclePlotter

warnings.simplefilter('ignore')

# ── cycle definition ───────────────────────────────────────────────────
WF          = 'R245fa'
P_COND      = 2.506e5        # Pa, condensing at 40 C
T_BRINE     = 403.15         # K, 130 C
M_BRINE     = 10.0           # kg/s
PPT         = 8.0            # K, specified minimum approach
DIV_N       = 200            # divisions used to build the profile
P_DESIGN    = 15e5           # Pa, the design point shown in panel (a)

# Pressures shown as markers in panel (b), in bar.
SWEEP_BAR   = [8, 10, 12, 15, 18, 20]

# The feasibility boundary in panel (b) is COMPUTED, not drawn by hand: it is
# the evaporating pressure at which the solver starts to refuse the
# specification (bisection between the last sampled pressure that solves and
# a higher one that is refused).
BOUNDARY_BRACKET = (SWEEP_BAR[-1], 21.0)   # bar; solves at the lower end
BOUNDARY_TOL     = 0.01                    # bar


def build(P_evap):
    """
    Build and solve the cycle at one evaporating pressure.

    The pump is solved first so that the evaporator's cold inlet enthalpy
    is known, leaving exactly two unknowns -- the brine outlet enthalpy
    and the working-fluid mass flow rate -- which the pinch condition and
    the energy balance close together.  It is solved again at the end,
    once the mass flow rate exists, so that its power can be evaluated.
    """
    M = ThermodynamicModel()
    M.set_dead_state()

    M.add_point(WF, '1', P=P_evap, Q=1)          # turbine inlet
    M.add_point(WF, '2', P=P_COND)               # turbine outlet
    M.add_point(WF, '3', P=P_COND, Q=0)          # condenser outlet
    M.add_point(WF, '4', P=P_evap)               # evaporator inlet
    M.add_point('water', 'b1', P=5e5, T=T_BRINE, Mass_flowrate=M_BRINE)
    M.add_point('water', 'b2', P=5e5)            # brine outlet unknown

    Pump(M, 'Pump', '3', '4', n_isen=0.75, Calculate=True)
    HeatExchanger(M, 'Evaporator', PPT=PPT, HEX_type='Evaporator',
                  HeatAdded=True,
                  Hot_In_state='b1', Hot_Out_state='b2',
                  Cold_In_state='4', Cold_Out_state='1', Calculate=True)
    Turbine(M, 'Turbine', '1', '2', n_isen=0.80, Calculate=True)
    HeatExchanger(M, 'Condenser', PPT=5, HEX_type='SimpleHEX',
                  HeatAdded=False,
                  Hot_In_state='2', Hot_Out_state='3',
                  Cold_In_state=None, Cold_Out_state=None, Calculate=True)
    Pump(M, 'Pump', '3', '4', n_isen=0.75, Calculate=True)

    M.ModelSummary(verbose=False)
    return M


def solves(p_bar):
    """True if the specification is accepted at this evaporating pressure."""
    try:
        build(p_bar * 1e5)
        return True
    except Exception:
        return False


def feasibility_limit(lo, hi, tol):
    """Highest pressure in [lo, hi] at which the solver still returns a result."""
    assert solves(lo) and not solves(hi), 'bracket must straddle the boundary'
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if solves(mid) else (lo, mid)
    return lo


# ── figure ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))

# (a) evaporator profile at the design point
M = build(P_DESIGN)
plotter = CyclePlotter(M)
min_dT = plotter.plot_hex_profile('Evaporator', div_N=DIV_N,
                                  save_csv=True, csv_path='orc_prof.csv')
prof = pd.read_csv('orc_prof.csv')
m_wf = M.Point['1'].Mass_flowrate

ax = axes[0]
Q = prof['Heat_Transferred_kW'] / 1e3
ax.plot(Q, prof['Hot_Temperature_C'],  '-', color='#c0392b', lw=2.0,
        label='geothermal brine')
ax.plot(Q, prof['Cold_Temperature_C'], '-', color='#1f6fb4', lw=2.0,
        label=WF)
i = prof['Delta_T_K'].idxmin()
ax.plot([Q[i], Q[i]],
        [prof['Cold_Temperature_C'][i], prof['Hot_Temperature_C'][i]],
        '--', color='#333', lw=1.1)
ax.annotate(f'pinch = {min_dT:.2f} K',
            xy=(Q[i], 0.5 * (prof['Hot_Temperature_C'][i]
                             + prof['Cold_Temperature_C'][i])),
            xytext=(8, -3), textcoords='offset points',
            fontsize=9.5, color='#2e7d32', fontweight='bold')
ax.set_xlabel('cumulative duty [MW]')
ax.set_ylabel('temperature [°C]')
ax.set_title(f'(a)  evaporator profile, {P_DESIGN/1e5:.0f} bar\n'
             f'$\\dot m$ = {m_wf:.3f} kg/s (derived)', fontsize=10)
ax.legend(frameon=False, fontsize=8.8, loc='lower right')
ax.grid(alpha=0.25)

# (b) efficiency and net power against evaporating pressure
rows = []
for p_bar in SWEEP_BAR:
    Mx = build(p_bar * 1e5)
    rows.append((p_bar, Mx.Efficiency, Mx.Net_power / 1e3,
                 Mx.Point['1'].Mass_flowrate,Mx.Point['b2'].T-273.15))
d = pd.DataFrame(rows, columns=['P', 'eta', 'W', 'm','T_brine_out_C'])

ax = axes[1]
ax.plot(d.P, d.eta, '-o', color='#2e7d32', lw=1.8, ms=6, label='efficiency')
ax.set_xlabel('evaporator pressure [bar]')
ax.set_ylabel('thermal efficiency [%]', color='#2e7d32')
ax.tick_params(axis='y', labelcolor='#2e7d32')

ax2 = ax.twinx()
ax2.plot(d.P, d.W, '-s', color='#1f6fb4', lw=1.8, ms=6, label='net power')
ax2.set_ylabel('net power [kW]', color='#1f6fb4')
ax2.tick_params(axis='y', labelcolor='#1f6fb4')

P_LIMIT = feasibility_limit(*BOUNDARY_BRACKET, BOUNDARY_TOL)
ax.axvline(P_LIMIT, color='#c0392b', ls=':', lw=1.6)
ax.text(P_LIMIT - 0.3, 0.52, f'refused above\n{P_LIMIT:.1f} bar\n(computed)', fontsize=8.2,
        color='#c0392b', ha='right', va='center', linespacing=1.3, transform=ax.get_xaxis_transform())

h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.8, loc='center left')
ax.set_title('(b)  design trade-off and its limit', fontsize=10)
ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig('fig_orc_example.png', dpi=300, bbox_inches='tight')

print(d.to_string(index=False))
print(f"computed feasibility limit at {PPT:g} K: {P_LIMIT:.2f} bar")
print(f"\ndesign point: m = {m_wf:.4f} kg/s, pinch = {min_dT:.2f} K ")
print("saved fig_orc_example.png")
