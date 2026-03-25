"""
main.py
-------
Regenerative Rankine cycle with open feedwater heater.
"""

from Thermosim import (
    ThermodynamicModel,
    Turbine, Pump, Splitter, Mixer, HeatExchanger,
)

# ── create model ──────────────────────────────────────────────
Model = ThermodynamicModel()
Model.set_dead_state()

wf    = 'water'
eta_t = 0.85
eta_p = 1.0
y     = 0.5          # initial guess for extraction fraction

# ── pressures ─────────────────────────────────────────────────
P1 = 8e6
P2 = 0.7e6
P3 = 0.008e6
P4 = P3
P5 = P2
P6 = P2
P7 = P1

# ── add state points ─────────────────────────────────────────
Model.add_point(wf, '1', P=P1, T=480+273.15, Mass_flowrate=1)
Model.add_point(wf, '2', P=P2)
Model.add_point(wf, '2a', P=P2)
Model.add_point(wf, '2b', P=P2)
Model.add_point(wf, '3', P=P3)
Model.add_point(wf, '4', P=P4, Q=0)
Model.add_point(wf, '5', P=P5)
Model.add_point(wf, '6', P=P6, Q=0)
Model.add_point(wf, '7', P=P7)

# ── first pass: get enthalpies we need ────────────────────────
Turbine(Model, ID="turbine1",
        In_state='1', Out_state='2',
        n_isen=eta_t, Calculate=True)

Splitter(Model, ID='Splitter',
         In_state='2', Out_states=['2a', '2b'],
         split_fractions=[y, 1-y], Calculate=True)

Pump(Model, ID='Pump1',
     In_state='4', Out_state='5',
     n_isen=eta_p, Calculate=True)

# ── compute correct extraction fraction ───────────────────────
y = (Model.Point['6'].H - Model.Point['5'].H) / \
    (Model.Point['2a'].H - Model.Point['5'].H)
print(f"Extraction fraction y = {y:.4f}")

# ── re-solve with correct y ──────────────────────────────────
Splitter(Model, ID='Splitter',
         In_state='2', Out_states=['2a', '2b'],
         split_fractions=[y, 1-y], Calculate=True)

Turbine(Model, ID="turbine2",
        In_state='2b', Out_state='3',
        n_isen=eta_t, Calculate=True)

Pump(Model, ID='Pump1',
     In_state='4', Out_state='5',
     n_isen=eta_p, Calculate=True)

HeatExchanger(Model, ID='Condenser',
              HEX_type='SimpleHEX', HeatAdded=False,
              Hot_In_state='3', Hot_Out_state='4',
              Cold_In_state=None, Cold_Out_state=None,
              PPT=5, Calculate=True)

Pump(Model, ID='Pump1',
     In_state='4', Out_state='5',
     n_isen=eta_p, Calculate=True)

Mixer(Model, ID='Mixer',
      In_states=['5', '2a'], Out_state='6',
      Calculate=True)

Pump(Model, ID='Pump2',
     In_state='6', Out_state='7',
     n_isen=eta_p, Calculate=True)

HeatExchanger(Model, ID='Boiler',
              PPT=5, HEX_type='SimpleHEX', HeatAdded=True,
              Hot_In_state=None, Hot_Out_state=None,
              Cold_In_state='7', Cold_Out_state='1',
              Calculate=True)

# ── print everything ─────────────────────────────────────────
print(Model)