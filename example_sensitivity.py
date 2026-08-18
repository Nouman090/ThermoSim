"""
example_sensitivity.py
----------------------
Demonstrates sensitivity analysis features.

Run with:   python example_sensitivity.py
"""

import numpy as np
from ThermoSim import (
    ThermodynamicModel, Turbine, Pump, HeatExchanger,
)
from ThermoSim.analysis import SensitivityAnalyzer


# ── Step 1: Define a BUILD FUNCTION ──────────────────────────
#
#    This function takes a dict of parameters, builds a complete
#    cycle, solves it, and returns the model.
#    The SensitivityAnalyzer will call this many times with
#    different parameter values.

def build_rankine(params):
    """Build and solve a simple Rankine cycle from parameters."""
    M = ThermodynamicModel()
    M.set_dead_state()

    T1    = params['T1']        # turbine inlet temp [K]
    P1    = params['P1']        # boiler pressure [Pa]
    P_cond = params['P_cond']   # condenser pressure [Pa]
    eta_t = params['eta_t']     # turbine isentropic efficiency
    eta_p = params['eta_p']     # pump isentropic efficiency

    M.add_point('water', '1', P=P1, T=T1, Mass_flowrate=1)
    M.add_point('water', '2', P=P_cond)
    M.add_point('water', '3', P=P_cond, Q=0)
    M.add_point('water', '4', P=P1)

    Turbine(M, 'Turbine', '1', '2', n_isen=eta_t, Calculate=True)
    HeatExchanger(M, 'Condenser', PPT=5, HEX_type='SimpleHEX',
                  HeatAdded=False,
                  Hot_In_state='2', Hot_Out_state='3',
                  Cold_In_state=None, Cold_Out_state=None,
                  Calculate=True)
    Pump(M, 'Pump', '3', '4', n_isen=eta_p, Calculate=True)
    HeatExchanger(M, 'Boiler', PPT=5, HEX_type='SimpleHEX',
                  HeatAdded=True,
                  Hot_In_state=None, Hot_Out_state=None,
                  Cold_In_state='4', Cold_Out_state='1',
                  Calculate=True)
    return M


# ── Step 2: Define base parameters ───────────────────────────
base = {
    'T1':     753.15,       # 480 °C
    'P1':     8e6,          # 80 bar
    'P_cond': 0.008e6,      # 0.08 bar
    'eta_t':  0.85,
    'eta_p':  1.0,
}

# ── Step 3: Create the analyzer ──────────────────────────────
sa = SensitivityAnalyzer(build_rankine, base)


# ── Example A: Single-parameter sweep ────────────────────────
print("=" * 60)
print("  SWEEP 1: Turbine inlet temperature")
print("=" * 60)

T_values = np.arange(623.15, 873.15, 25)    # 350–600 °C

df1 = sa.single_sweep(
    param_name='T1',
    values=T_values,
    outputs={
        'Efficiency (%)':     lambda m: m.Efficiency,
        'Net Power (kW)':     lambda m: m.Net_power / 1e3,
    },
    x_label='Turbine Inlet Temperature [K]',
)

print("\nResults table:")
print(df1.to_string(index=False))


# ── Example B: Sweep turbine efficiency ──────────────────────
print("\n" + "=" * 60)
print("  SWEEP 2: Turbine isentropic efficiency")
print("=" * 60)

eta_values = np.arange(0.70, 1.01, 0.05)

df2 = sa.single_sweep(
    param_name='eta_t',
    values=eta_values,
    outputs={
        'Efficiency (%)': lambda m: m.Efficiency,
    },
    x_label='Turbine Isentropic Efficiency',
)


# ── Example C: Multi-output on one chart ─────────────────────
print("\n" + "=" * 60)
print("  SWEEP 3: Multi-output vs boiler pressure")
print("=" * 60)

P_values = np.arange(4e6, 16e6, 1e6)

df3 = sa.multi_output_sweep(
    param_name='P1',
    values=P_values,
    outputs={
        'Efficiency (%)':     lambda m: m.Efficiency,
        'Net Power (kW)':     lambda m: m.Net_power / 1e3,
        'Turbine Work (kW)':  lambda m: m.Component['Turbine'].work / 1e3,
    },
    x_label='Boiler Pressure [Pa]',
)


# ── Example D: Two-parameter contour ─────────────────────────
print("\n" + "=" * 60)
print("  SWEEP 4: Contour — T1 vs P1 → Efficiency")
print("=" * 60)

T_vals = np.arange(623.15, 873.15, 50)
P_vals = np.arange(4e6, 14e6, 2e6)

df4 = sa.double_sweep(
    param1_name='T1', param1_values=T_vals,
    param2_name='P1', param2_values=P_vals,
    output_name='Efficiency (%)',
    output_func=lambda m: m.Efficiency,
)

print("\nDone! Close plot windows to exit.")