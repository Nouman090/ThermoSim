"""
test_heat_exchanger_2unknown_bracket.py
-----------------------------------------
Regression test for a bracket bug in HeatExchanger._solve_2_unknown's
{H_hi,mh} / {H_ho,mh} / {H_ci,mc} cases.

Each of these cases solves for a free enthalpy H_x via
mass_ = Q/(H_x - H_fixed), which has a POLE at H_x == H_fixed: on one side
of it mass_ is physically valid (> 0), on the other it is negative (the
"hot" or "cold" stream running backwards). The original bracket spanned
straight across that pole (e.g. H_lo = a generic lower bound well below
H_fixed), so brentq could converge on the unphysical side -- silently
returning a NEGATIVE Hot_Mass_flowrate/Cold_Mass_flowrate that still zeroes
the residual arithmetically, with nothing downstream flagging it as wrong.

The fix tightens the bracket to stay strictly on the physical side of the
pole. These tests reproduce the original failure (negative mass flow, hot
inlet colder than hot outlet) and confirm it's gone, for both a plain
CoolProp fluid and an NH3-H2O mixture.

Run with:   pytest tests/test_heat_exchanger_2unknown_bracket.py -v
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


class TestHhiMhBracket:
    """{'H_hi', 'mh'}: pole at H_hi == Hot_Out.H."""

    def test_pure_fluid_finds_physical_root(self):
        M = ThermodynamicModel()
        M.add_point('Air', 'hi', P=1e5)             # H, mh unknown
        M.add_point('Air', 'ho', P=1e5, T=350.0)     # fixed
        M.add_point('Nitrogen', 'ci', P=1e5, T=300.0, Mass_flowrate=2.0)
        M.add_point('Nitrogen', 'co', P=1e5, T=330.0, Mass_flowrate=2.0)
        hx = HeatExchanger(M, 'HX1', PPT=30.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Hot_Mass_flowrate > 0
        assert hx.Hot_In.T > hx.Hot_Out.T
        assert hx.min_dT == pytest.approx(hx.PPT, abs=1e-1)

    def test_mixture_finds_physical_root(self):
        M = ThermodynamicModel()
        M.add_point(NH3H2OFluid(0.995), 'hi', P=15e5)
        M.add_point(NH3H2OFluid(0.995), 'ho', P=15e5, Q=0.0)
        M.add_point('water', 'ci', P=2e5, T=288.15, Mass_flowrate=1.0)
        M.add_point('water', 'co', P=2e5, T=300.0, Mass_flowrate=1.0)
        hx = HeatExchanger(M, 'CondB', PPT=15.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Hot_Mass_flowrate > 0
        assert hx.Hot_In.T > hx.Hot_Out.T
        assert hx.min_dT == pytest.approx(hx.PPT, abs=1e-1)

    def test_infeasible_target_raises_instead_of_wrong_answer(self):
        """A PPT tighter than the branch's achievable minimum must fail
        loudly, not return a negative-mass-flow 'solution'."""
        M = ThermodynamicModel()
        M.add_point('Ammonia', 'hi', P=15e5)
        M.add_point('Ammonia', 'ho', P=15e5, Q=0.0)
        M.add_point('water', 'wi', P=2e5, T=298.15, Mass_flowrate=11.0)
        M.add_point('water', 'wo', P=2e5, T=305.0, Mass_flowrate=11.0)
        with pytest.raises(Exception):
            HeatExchanger(M, 'CondC', PPT=5.0, HEX_type='Condenser',
                          HeatAdded=False,
                          Hot_In_state='hi', Hot_Out_state='ho',
                          Cold_In_state='wi', Cold_Out_state='wo',
                          Calculate=True)


class TestHhoMhBracket:
    """{'H_ho', 'mh'}: pole at H_ho == Hot_In.H (mirror case)."""

    def test_pure_fluid_finds_physical_root(self):
        M = ThermodynamicModel()
        M.add_point('Ammonia', 'hi', P=15e5, Q=1.0)   # fixed
        M.add_point('Ammonia', 'ho', P=15e5)           # H, mh unknown
        M.add_point('water', 'wi', P=2e5, T=288.15, Mass_flowrate=1.0)
        M.add_point('water', 'wo', P=2e5, T=300.0, Mass_flowrate=1.0)
        hx = HeatExchanger(M, 'CondC', PPT=5.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='wi', Cold_Out_state='wo',
                           Calculate=True)
        assert hx.Hot_Mass_flowrate > 0
        assert hx.Hot_Out.T < hx.Hot_In.T


class TestHciMcBracket:
    """{'H_ci', 'mc'}: pole at H_ci == Cold_Out.H."""

    def test_pure_fluid_finds_physical_root(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, T=350.0, Mass_flowrate=1.0)
        M.add_point('water', 'ho', P=2e5, T=310.0, Mass_flowrate=1.0)
        M.add_point('Nitrogen', 'ci', P=1e5)            # H, mc unknown
        M.add_point('Nitrogen', 'co', P=1e5, T=300.0)   # fixed
        hx = HeatExchanger(M, 'HX2', PPT=15.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Cold_Mass_flowrate > 0
        assert hx.Cold_In.T < hx.Cold_Out.T
        assert hx.min_dT == pytest.approx(hx.PPT, abs=1e-1)
