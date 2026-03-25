"""
test_integration.py
-------------------
End-to-end tests: build a complete cycle, solve, check results.

Run with:   pytest tests/test_integration.py -v
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Thermosim import (
    ThermodynamicModel,
    Turbine, Pump, Compressor,
    Splitter, Mixer,
    HeatExchanger,
)


def fresh_model():
    m = ThermodynamicModel()
    m.set_dead_state(T0=298.15, P0=101325)
    return m


# ================================================================== #
#  Simple Rankine Cycle
# ================================================================== #
class TestSimpleRankine:
    """
    1 → Turbine → 2 → Condenser → 3 → Pump → 4 → Boiler → 1
    """

    def _build(self, eta_t=0.85, eta_p=1.0, T1=753.15, P1=8e6, P_cond=0.008e6):
        M = fresh_model()
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

    def test_all_components_solved(self):
        M = self._build()
        for comp in M.Component.values():
            assert comp.Solution_Status is True, f"{comp.ID} not solved"

    def test_positive_efficiency(self):
        M = self._build()
        M.ModelSummary()
        assert M.Efficiency > 0
        assert M.Efficiency < 100

    def test_net_power_positive(self):
        M = self._build()
        M.ModelSummary()
        assert M.Net_power > 0

    def test_energy_balance(self):
        """Q_in = W_net + Q_out (first law)."""
        M = self._build()
        M.ModelSummary()
        # Q_out is heat rejected (positive number)
        # Energy balance: Q_in = W_net + Q_out
        balance = M.Q_in - M.Net_power - M.Q_out
        assert abs(balance) < 1.0   # within 1 W

    def test_higher_temp_better_efficiency(self):
        """Higher turbine inlet T → higher efficiency."""
        M1 = self._build(T1=673.15)
        M1.ModelSummary()
        M2 = self._build(T1=873.15)
        M2.ModelSummary()
        assert M2.Efficiency > M1.Efficiency

    def test_isentropic_better_than_real(self):
        """Isentropic turbine gives better efficiency."""
        M1 = self._build(eta_t=1.0)
        M1.ModelSummary()
        M2 = self._build(eta_t=0.80)
        M2.ModelSummary()
        assert M1.Efficiency > M2.Efficiency


# ================================================================== #
#  Regenerative Rankine (with extraction)
# ================================================================== #
class TestRegenerativeRankine:
    """
    The cycle from the user's example:
    1 → Turbine1 → 2 → Splitter → 2a (extraction) + 2b
    2b → Turbine2 → 3 → Condenser → 4 → Pump1 → 5
    Mixer(5 + 2a) → 6 → Pump2 → 7 → Boiler → 1
    """

    def _build(self):
        M = fresh_model()
        wf = 'water'
        eta_t = 0.85
        eta_p = 1.0

        P1 = 8e6; P2 = 0.7e6; P3 = 0.008e6
        P4 = P3;  P5 = P2;    P6 = P2; P7 = P1

        M.add_point(wf, '1', P=P1, T=480+273.15, Mass_flowrate=1)
        M.add_point(wf, '2', P=P2)
        M.add_point(wf, '2a', P=P2)
        M.add_point(wf, '2b', P=P2)
        M.add_point(wf, '3', P=P3)
        M.add_point(wf, '4', P=P4, Q=0)
        M.add_point(wf, '5', P=P5)
        M.add_point(wf, '6', P=P6, Q=0)
        M.add_point(wf, '7', P=P7)

        # first pass
        Turbine(M, 'turbine1', '1', '2', n_isen=eta_t, Calculate=True)
        Splitter(M, 'Splitter', '2', ['2a', '2b'], [0.5, 0.5], Calculate=True)
        Pump(M, 'Pump1', '4', '5', n_isen=eta_p, Calculate=True)

        # correct extraction fraction
        y = ((M.Point['6'].H - M.Point['5'].H) /
             (M.Point['2a'].H - M.Point['5'].H))

        Splitter(M, 'Splitter', '2', ['2a', '2b'], [y, 1-y], Calculate=True)
        Turbine(M, 'turbine2', '2b', '3', n_isen=eta_t, Calculate=True)

        HeatExchanger(M, 'Condenser', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=False,
                      Hot_In_state='3', Hot_Out_state='4',
                      Cold_In_state=None, Cold_Out_state=None,
                      Calculate=True)

        Pump(M, 'Pump1', '4', '5', n_isen=eta_p, Calculate=True)
        Mixer(M, 'Mixer', ['5', '2a'], '6', Calculate=True)
        Pump(M, 'Pump2', '6', '7', n_isen=eta_p, Calculate=True)

        HeatExchanger(M, 'Boiler', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=True,
                      Hot_In_state=None, Hot_Out_state=None,
                      Cold_In_state='7', Cold_Out_state='1',
                      Calculate=True)
        return M

    def test_regenerative_solves(self):
        M = self._build()
        # all key components solved
        for key in ['turbine1', 'turbine2', 'Pump1', 'Pump2',
                     'Condenser', 'Boiler', 'Mixer']:
            assert M.Component[key].Solution_Status is True

    def test_regenerative_efficiency(self):
        M = self._build()
        M.ModelSummary()
        # regenerative should be > 25% efficient
        assert M.Efficiency > 25

    def test_mass_conservation_at_mixer(self):
        M = self._build()
        m_2a = M.Point['2a'].Mass_flowrate
        m_5  = M.Point['5'].Mass_flowrate
        m_6  = M.Point['6'].Mass_flowrate
        assert abs(m_2a + m_5 - m_6) < 0.001

    def test_mass_conservation_at_splitter(self):
        M = self._build()
        m_2  = M.Point['2'].Mass_flowrate
        m_2a = M.Point['2a'].Mass_flowrate
        m_2b = M.Point['2b'].Mass_flowrate
        assert abs(m_2 - m_2a - m_2b) < 0.001


# ================================================================== #
#  Simple Brayton Cycle  (uses Compressor)
# ================================================================== #
class TestSimpleBrayton:
    """
    1 → Compressor → 2 → Heater → 3 → Turbine → 4 → Cooler → 1
    """

    def _build(self, T1=300, P1=1e5, rp=10, T3=1200, eta_c=0.85, eta_t=0.85):
        M = fresh_model()
        P2 = P1 * rp

        M.add_point('Air', '1', P=P1, T=T1, Mass_flowrate=1)
        M.add_point('Air', '2', P=P2)
        M.add_point('Air', '3', P=P2, T=T3)
        M.add_point('Air', '4', P=P1)

        Compressor(M, 'Comp', '1', '2', n_isen=eta_c, Calculate=True)
        HeatExchanger(M, 'Heater', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=True,
                      Hot_In_state=None, Hot_Out_state=None,
                      Cold_In_state='2', Cold_Out_state='3',
                      Calculate=True)
        Turbine(M, 'Turb', '3', '4', n_isen=eta_t, Calculate=True)
        HeatExchanger(M, 'Cooler', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=False,
                      Hot_In_state='4', Hot_Out_state='1',
                      Cold_In_state=None, Cold_Out_state=None,
                      Calculate=True)
        return M

    def test_brayton_solves(self):
        M = self._build()
        for comp in M.Component.values():
            assert comp.Solution_Status is True

    def test_brayton_positive_efficiency(self):
        M = self._build()
        M.ModelSummary()
        assert M.Efficiency > 0

    def test_turbine_produces_more_than_compressor(self):
        M = self._build()
        M.ModelSummary()
        assert M.Net_power > 0

    def test_compressor_pressure_ratio(self):
        M = self._build(rp=12)
        assert abs(M.Component['Comp'].pressure_ratio - 12.0) < 0.1