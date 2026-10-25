"""
test_heat_exchanger_same_side_pair.py
--------------------------------------
Regression tests for {H_hi, H_ho} (case 1) and {H_ci, H_co} (case 2) in
HeatExchanger._solve_2_unknown.

These were previously rejected outright as "Q undefined on that side", but
whenever exactly these two are unknown the OTHER side is fully known by
definition, so Q is already fixed from it. The unknown corner is searched
directly -- it IS the 'cold_to_hot' marching start point for its stream --
with no assumption about where the pinch sits, unlike an anchor. That
matters for a phase-changing stream, whose true pinch can sit inside the
unit rather than at either end.

Run with:   pytest tests/test_heat_exchanger_same_side_pair.py -v
"""
import sys
import os
import warnings

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, NH3H2OFluid, HeatExchanger


@pytest.fixture(autouse=True)
def _quiet():
    warnings.filterwarnings('ignore')


def _pinch_tol(PPT):
    return max(1e-2, 1e-3 * PPT)


class TestHhiHhoBothHotUnknown:
    """{'H_hi', 'H_ho'} (case 1): cold side + mh fully known."""

    def test_phase_changing_hot_stream(self):
        """Hot stream condenses across the unit, so the true pinch sits at
        an internal saturation crossing -- exactly where an anchor-based
        (pinch-at-an-end) approach would get the wrong answer."""
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, Mass_flowrate=2.0)   # H_hi, H_ho unknown
        M.add_point('water', 'ho', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=13.553)
        M.add_point('water', 'co', P=3e5, T=387.35, Mass_flowrate=13.553)  # ~114.2 C
        hx = HeatExchanger(M, 'HX1', PPT=8.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
        assert hx.Hot_In.T > hx.Hot_Out.T > M.Point['ci'].T
        Q_hot = hx.Hot_Mass_flowrate * (hx.Hot_In.H - hx.Hot_Out.H)
        Q_cold = 13.553 * (M.Point['co'].H - M.Point['ci'].H)
        assert Q_hot == pytest.approx(Q_cold, rel=1e-3)

    def test_single_phase_hot_stream(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ho', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=10.0)
        M.add_point('water', 'co', P=3e5, T=340.0, Mass_flowrate=10.0)
        hx = HeatExchanger(M, 'HX2', PPT=8.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
        assert hx.Hot_In.T > hx.Hot_Out.T

    def test_condenser_type(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ho', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=10.0)
        M.add_point('water', 'co', P=3e5, T=340.0, Mass_flowrate=10.0)
        hx = HeatExchanger(M, 'HX3', PPT=8.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))

    def test_evaporator_type(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ho', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=10.0)
        M.add_point('water', 'co', P=3e5, T=340.0, Mass_flowrate=10.0)
        hx = HeatExchanger(M, 'HX4', PPT=8.0, HEX_type='Evaporator',
                           HeatAdded=True,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))

    def test_different_fluids_each_side(self):
        M = ThermodynamicModel()
        M.add_point('Air', 'hi', P=2e5, Mass_flowrate=2.0)
        M.add_point('Air', 'ho', P=2e5, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=0.3)
        M.add_point('water', 'co', P=3e5, T=340.0, Mass_flowrate=0.3)
        hx = HeatExchanger(M, 'HX5', PPT=8.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))

    def test_tight_ppt(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, Mass_flowrate=1.5)
        M.add_point('water', 'ho', P=2e5, Mass_flowrate=1.5)
        M.add_point('water', 'ci', P=3e5, T=320.0, Mass_flowrate=6.25)
        M.add_point('water', 'co', P=3e5, T=360.0, Mass_flowrate=6.25)
        hx = HeatExchanger(M, 'HX6', PPT=0.5, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))

    def test_nh3h2o_mixture(self):
        M = ThermodynamicModel()
        M.add_point(NH3H2OFluid(0.995), 'hi', P=15e5, Mass_flowrate=1.0)
        M.add_point(NH3H2OFluid(0.995), 'ho', P=15e5, Mass_flowrate=1.0)
        M.add_point('water', 'ci', P=2e5, T=288.15, Mass_flowrate=1.0)
        M.add_point('water', 'co', P=2e5, T=300.0, Mass_flowrate=1.0)
        hx = HeatExchanger(M, 'CondM', PPT=5.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=1e-1)


class TestHciHcoBothColdUnknown:
    """{'H_ci', 'H_co'} (case 2): hot side + mc fully known (mirror case)."""

    def test_phase_changing_hot_stream(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, T=420.0, Mass_flowrate=2.0)
        M.add_point('water', 'ho', P=2e5, T=340.0, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, Mass_flowrate=13.553)   # H_ci, H_co unknown
        M.add_point('water', 'co', P=3e5, Mass_flowrate=13.553)
        hx = HeatExchanger(M, 'HX7', PPT=8.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
        assert hx.Cold_In.T < hx.Cold_Out.T
        Q_hot = 2.0 * (M.Point['hi'].H - M.Point['ho'].H)
        Q_cold = hx.Cold_Mass_flowrate * (hx.Cold_Out.H - hx.Cold_In.H)
        assert Q_hot == pytest.approx(Q_cold, rel=1e-3)

    def test_single_phase(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, T=380.0, Mass_flowrate=2.0)
        M.add_point('water', 'ho', P=2e5, T=350.0, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, Mass_flowrate=10.0)
        M.add_point('water', 'co', P=3e5, Mass_flowrate=10.0)
        hx = HeatExchanger(M, 'HX8', PPT=8.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))

    def test_different_fluids_each_side(self):
        M = ThermodynamicModel()
        M.add_point('Air', 'hi', P=2e5, T=420.0, Mass_flowrate=2.0)
        M.add_point('Air', 'ho', P=2e5, T=340.0, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, Mass_flowrate=0.3)
        M.add_point('water', 'co', P=3e5, Mass_flowrate=0.3)
        hx = HeatExchanger(M, 'HX9', PPT=8.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))


class TestStillUnsolvableCasesUnaffected:
    """Cases 8, 9 and 15 must still raise -- only 1 and 2 changed."""

    def test_case_8_still_raises(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, Mass_flowrate=2.0)        # H_hi unknown
        M.add_point('water', 'ho', P=2e5, T=340.0, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=10.0)
        M.add_point('water', 'co', P=3e5, Mass_flowrate=10.0)       # H_co unknown
        with pytest.raises(ValueError, match="case 8"):
            HeatExchanger(M, 'HX10', PPT=8.0, HEX_type='double_pipe',
                          HeatAdded=None,
                          Hot_In_state='hi', Hot_Out_state='ho',
                          Cold_In_state='ci', Cold_Out_state='co',
                          Calculate=True)

    def test_case_9_still_raises(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, T=420.0, Mass_flowrate=2.0)
        M.add_point('water', 'ho', P=2e5, Mass_flowrate=2.0)        # H_ho unknown
        M.add_point('water', 'ci', P=3e5, Mass_flowrate=10.0)       # H_ci unknown
        M.add_point('water', 'co', P=3e5, T=390.0, Mass_flowrate=10.0)
        with pytest.raises(ValueError, match="case 9"):
            HeatExchanger(M, 'HX11', PPT=8.0, HEX_type='double_pipe',
                          HeatAdded=None,
                          Hot_In_state='hi', Hot_Out_state='ho',
                          Cold_In_state='ci', Cold_Out_state='co',
                          Calculate=True)

    def test_case_15_still_raises(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, T=420.0)
        M.add_point('water', 'ho', P=2e5, T=340.0)
        M.add_point('water', 'ci', P=3e5, T=300.0)
        M.add_point('water', 'co', P=3e5, T=390.0)
        with pytest.raises(ValueError, match="case 15"):
            HeatExchanger(M, 'HX12', PPT=8.0, HEX_type='double_pipe',
                          HeatAdded=None,
                          Hot_In_state='hi', Hot_Out_state='ho',
                          Cold_In_state='ci', Cold_Out_state='co',
                          Calculate=True)
