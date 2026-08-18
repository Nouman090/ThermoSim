"""
ThermoSim verification against published textbook worked examples.

Reference
---------
Moran, M. J. and Shapiro, H. N., *Fundamentals of Engineering
Thermodynamics*, 5th edition (SI units).  Examples 8.1 - 8.4.

Method
------
Each cycle is built in ThermoSim exactly as specified in the source, solved
on a unit mass-flow basis, and the mass flow is then scaled to the net power
the textbook specifies.  Published and computed values are compared as a
relative error.

Note on property data: Moran's tables predate IAPWS-95, which is what
CoolProp implements.  Agreement of the raw property values is ~0.04 %, so
that is the floor on achievable agreement for any derived quantity.
"""
import math
from ThermoSim import ThermodynamicModel, Turbine, Pump, HeatExchanger

W_NET_TARGET = 100e6          # W, common to all four examples


# --------------------------------------------------------------------- #
#  Cycle builders  (unit mass flow; scaled afterwards)
# --------------------------------------------------------------------- #
def rankine(eta_t=1.0, eta_p=1.0):
    """Examples 8.1 and 8.2: saturated vapour at 8 MPa, condensing at 8 kPa."""
    M = ThermodynamicModel()
    M.set_dead_state()
    M.add_point('water', '1', P=8e6,     Q=1, Mass_flowrate=1)
    M.add_point('water', '2', P=0.008e6)
    M.add_point('water', '3', P=0.008e6, Q=0)
    M.add_point('water', '4', P=8e6)

    Turbine(M, 'Turbine', '1', '2', n_isen=eta_t, Calculate=True)
    HeatExchanger(M, 'Condenser', PPT=5, HEX_type='SimpleHEX', HeatAdded=False,
                  Hot_In_state='2', Hot_Out_state='3',
                  Cold_In_state=None, Cold_Out_state=None, Calculate=True)
    Pump(M, 'Pump', '3', '4', n_isen=eta_p, Calculate=True)
    HeatExchanger(M, 'Boiler', PPT=5, HEX_type='SimpleHEX', HeatAdded=True,
                  Hot_In_state=None, Hot_Out_state=None,
                  Cold_In_state='4', Cold_Out_state='1', Calculate=True)
    return M, ['1', '2', '3', '4']


def reheat(eta_t=1.0, eta_p=1.0):
    """Examples 8.3 and 8.4: 8 MPa / 480 C, reheat at 0.7 MPa to 440 C."""
    M = ThermodynamicModel()
    M.set_dead_state()
    M.add_point('water', '1', P=8e6,     T=753.15, Mass_flowrate=1)   # 480 C
    M.add_point('water', '2', P=0.7e6)
    M.add_point('water', '3', P=0.7e6,   T=713.15)                    # 440 C
    M.add_point('water', '4', P=0.008e6)
    M.add_point('water', '5', P=0.008e6, Q=0)
    M.add_point('water', '6', P=8e6)

    Turbine(M, 'HP Turbine', '1', '2', n_isen=eta_t, Calculate=True)
    HeatExchanger(M, 'Reheater', PPT=5, HEX_type='SimpleHEX', HeatAdded=True,
                  Hot_In_state=None, Hot_Out_state=None,
                  Cold_In_state='2', Cold_Out_state='3', Calculate=True)
    Turbine(M, 'LP Turbine', '3', '4', n_isen=eta_t, Calculate=True)
    HeatExchanger(M, 'Condenser', PPT=5, HEX_type='SimpleHEX', HeatAdded=False,
                  Hot_In_state='4', Hot_Out_state='5',
                  Cold_In_state=None, Cold_Out_state=None, Calculate=True)
    Pump(M, 'Pump', '5', '6', n_isen=eta_p, Calculate=True)
    HeatExchanger(M, 'Boiler', PPT=5, HEX_type='SimpleHEX', HeatAdded=True,
                  Hot_In_state=None, Hot_Out_state=None,
                  Cold_In_state='6', Cold_Out_state='1', Calculate=True)
    return M, ['1', '2', '3', '4', '5', '6']


# --------------------------------------------------------------------- #
#  Reference data, transcribed from the source
# --------------------------------------------------------------------- #
CASES = [
    dict(no='8.1', name='Ideal Rankine',
         build=lambda: rankine(1.0, 1.0),
         eta=37.1, m_dot_kgph=3.77e5,
         h=[2758.0, 1794.8, 173.88, 181.94]),
    dict(no='8.2', name='Rankine with irreversibilities',
         build=lambda: rankine(0.85, 0.85),
         eta=31.4, m_dot_kgph=4.449e5,
         h=[2758.0, 1939.3, 173.88, 183.36]),
    dict(no='8.3', name='Ideal reheat',
         build=lambda: reheat(1.0, 1.0),
         eta=40.3, m_dot_kgph=2.363e5,
         h=[3348.4, 2741.8, 3353.3, 2428.5, 173.88, 181.94]),
    dict(no='8.4', name='Reheat, turbine irreversibility',
         build=lambda: reheat(0.85, 1.0),
         eta=35.1, m_dot_kgph=None,
         h=[3348.4, 2832.8, 3353.3, 2567.2, 173.88, 181.94]),
]


def err(computed, published):
    return 100.0 * (computed - published) / published


def run():
    print("=" * 78)
    print("ThermoSim verification -- Moran & Shapiro, 5th ed. (SI), Ex. 8.1-8.4")
    print("=" * 78)

    summary = []
    for c in CASES:
        M, names = c['build']()
        M.ModelSummary(verbose=False)

        # scale unit-mass results to the textbook's 100 MW net output
        w_net_specific = M.Net_power                 # W per (kg/s)
        m_dot = W_NET_TARGET / w_net_specific        # kg/s
        m_dot_kgph = m_dot * 3600.0

        print(f"\n--- Example {c['no']}: {c['name']} ---")
        print(f"  {'':22s} {'published':>12s} {'ThermoSim':>12s} {'error %':>9s}")
        print(f"  {'thermal efficiency %':22s} {c['eta']:12.2f} "
              f"{M.Efficiency:12.2f} {err(M.Efficiency, c['eta']):9.2f}")
        if c['m_dot_kgph'] is not None:
            print(f"  {'mass flow kg/h':22s} {c['m_dot_kgph']:12.4g} "
                  f"{m_dot_kgph:12.4g} {err(m_dot_kgph, c['m_dot_kgph']):9.2f}")

        worst_h = 0.0
        for i, (name, h_pub) in enumerate(zip(names, c['h']), start=1):
            h_com = M.Point[name].H / 1e3
            e = err(h_com, h_pub)
            worst_h = max(worst_h, abs(e))
            print(f"  {'h' + str(i) + ' kJ/kg':22s} {h_pub:12.2f} "
                  f"{h_com:12.2f} {e:9.2f}")

        summary.append((c['no'], c['name'], c['eta'], M.Efficiency,
                        err(M.Efficiency, c['eta']), worst_h))

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"{'Ex.':5s} {'cycle':32s} {'pub %':>7s} {'ThermoSim %':>12s} "
          f"{'err %':>7s} {'max |h| err %':>14s}")
    for no, name, pub, com, e, wh in summary:
        print(f"{no:5s} {name:32s} {pub:7.1f} {com:12.2f} {e:7.2f} {wh:14.2f}")

    worst = max(abs(s[4]) for s in summary)
    print(f"\nLargest efficiency deviation: {worst:.2f} %")


if __name__ == '__main__':
    run()
