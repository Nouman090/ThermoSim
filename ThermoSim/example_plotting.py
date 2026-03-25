"""
example_plotting.py
-------------------
Demonstrates all plotting features.

Run with:   python example_plotting.py
"""

from Thermosim import (
    ThermodynamicModel, Turbine, Pump, HeatExchanger,
)
from Thermosim.plotting import CyclePlotter

# ── Build a simple Rankine cycle ──────────────────────────────
M = ThermodynamicModel()
M.set_dead_state()

M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
M.add_point('water', '2', P=0.008e6)
M.add_point('water', '3', P=0.008e6, Q=0)
M.add_point('water', '4', P=8e6)

Turbine(M, 'Turbine', '1', '2', n_isen=0.85, Calculate=True)

HeatExchanger(M, 'Condenser', PPT=5, HEX_type='SimpleHEX',
              HeatAdded=False,
              Hot_In_state='2', Hot_Out_state='3',
              Cold_In_state=None, Cold_Out_state=None,
              Calculate=True)

Pump(M, 'Pump', '3', '4', n_isen=1.0, Calculate=True)

HeatExchanger(M, 'Boiler', PPT=5, HEX_type='SimpleHEX',
              HeatAdded=True,
              Hot_In_state=None, Hot_Out_state=None,
              Cold_In_state='4', Cold_Out_state='1',
              Calculate=True)

M.ModelSummary()

# ── create the plotter ────────────────────────────────────────
plotter = CyclePlotter(M)

# 1. T-s Diagram
print("\n--- T-s Diagram ---")
plotter.plot_Ts_diagram(['1', '2', '3', '4', '1'])

# 2. P-h Diagram
print("\n--- P-h Diagram ---")
plotter.plot_Ph_diagram(['1', '2', '3', '4', '1'])

# 3. h-s Diagram
print("\n--- h-s Diagram ---")
plotter.plot_hs_diagram(['1', '2', '3', '4', '1'])

# 4. Exergy bar chart
print("\n--- Exergy Destruction Bar Chart ---")
plotter.plot_exergy_bar()

# 5. Exergy pie chart
print("\n--- Exergy Destruction Pie Chart ---")
plotter.plot_exergy_pie()

# 6. Energy summary
print("\n--- Energy Summary ---")
plotter.plot_energy_summary()

print("\nDone! Close plot windows to exit.")