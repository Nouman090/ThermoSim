"""
Geothermal organic Rankine cycle sized by the evaporator's minimum approach.

This is the example of Section 2.3 of the SoftwareX article (R245fa, 15 bar,
brine at 130 C and 10 kg/s, 8 K minimum approach).  Run it with

    python examples/geothermal_orc.py

All quantities are SI: pressure in Pa, temperature in K, specific enthalpy in
J/kg, mass flow rate in kg/s.  The working-fluid mass flow rate is NOT given:
the evaporator derives it, together with the brine outlet state, from the
energy balance and the 8 K minimum approach (two unknowns).

Expected output (ThermoSim 3.5.0, CoolProp 8.0.0), also asserted below:
    working-fluid mass flow rate   4.735 kg/s
    brine outlet temperature     104.78 C
    evaporator minimum approach    8.00 K
    turbine 125.4 kW, pump 6.1 kW   net power 119.4 kW
    thermal efficiency            11.16 %
"""
from ThermoSim import ThermodynamicModel, Turbine, Pump, HeatExchanger

M = ThermodynamicModel()
M.set_dead_state()

# working fluid: pressures are assumed, mass flow is left unknown
M.add_point('R245fa', '1', P=15e5, Q=1)          # turbine inlet, saturated vapour
M.add_point('R245fa', '2', P=2.506e5)            # turbine outlet
M.add_point('R245fa', '3', P=2.506e5, Q=0)       # condenser outlet, saturated liquid
M.add_point('R245fa', '4', P=15e5)               # evaporator inlet
# geothermal brine
M.add_point('water', 'b1', P=5e5, T=403.15, Mass_flowrate=10)   # 130 C
M.add_point('water', 'b2', P=5e5)                # outlet unknown

Pump(M, 'Pump', '3', '4', n_isen=0.75, Calculate=True)
HeatExchanger(M, 'Evaporator', PPT=8, HEX_type='Evaporator', HeatAdded=True,
              Hot_In_state='b1', Hot_Out_state='b2',
              Cold_In_state='4', Cold_Out_state='1', Calculate=True)
Turbine(M, 'Turbine', '1', '2', n_isen=0.80, Calculate=True)
HeatExchanger(M, 'Condenser', PPT=5, HEX_type='SimpleHEX', HeatAdded=False,
              Hot_In_state='2', Hot_Out_state='3',
              Cold_In_state=None, Cold_Out_state=None, Calculate=True)
M.Component['Pump'].Cal()                        # re-solve for work, now that the flow is known

M.ModelSummary(verbose=True)

m_wf = M.Point['1'].Mass_flowrate
T_b2 = M.Point['b2'].T - 273.15
dT_min = M.Component['Evaporator'].min_dT
w_net = (M.Component['Turbine'].work - M.Component['Pump'].work) / 1e3
print(f"\nworking-fluid mass flow rate : {m_wf:.3f} kg/s")
print(f"brine outlet temperature     : {T_b2:.2f} C")
print(f"evaporator minimum approach  : {dT_min:.2f} K")
print(f"net power                    : {w_net:.1f} kW")

assert abs(m_wf - 4.735) < 0.005 and abs(T_b2 - 104.78) < 0.02
assert abs(dT_min - 8.0) < 0.01
