"""
test_heat_exchanger_3unknown_anchor_sanity.py
-----------------------------------------------
Regression test for a bug in HeatExchanger._apply_anchor's ordering
sanity checks, found via an independent benchmark (anchor_sanity_check.patch).

Two of the four PPT-diagonal anchors sanity-check the resulting enthalpy
against the OPPOSITE corner's temperature before accepting it:
  - H_ho -> H_ci (Condenser) requires Cold_Out.T
  - H_co -> H_hi (Evaporator) requires Hot_Out.T / Hot_In.T respectively

For a 3-unknown reduction, that opposite corner can itself be one of the
three unknowns -- it isn't resolved until the *remaining* 2-unknown solve
runs, afterward. Requiring it up front made every such combination raise
"temperature is unknown", even though a valid anchor existed and the
remaining 2-unknown case was already solvable on its own:

  Evaporator: {H_hi, H_ho, mh} and {H_hi, H_ho, mc}
  Condenser:  {H_ci, H_co, mh} and {H_ci, H_co, mc}

double_pipe was never affected: for these same combinations it has an
unaffected alternate anchor (H_ci->H_ho / H_hi->H_co) that the anchor
list tries first and which carries no such check, so it always
succeeded. Evaporator and Condenser each have only one of the two
diagonals available, so when that one was blocked, they had nothing to
fall back to.

Fixed by only enforcing the ordering check when the opposite corner's
temperature is already known; otherwise defer to the downstream pinch
and energy-balance checks that already run on every solve.

Run with:   pytest tests/test_heat_exchanger_3unknown_anchor_sanity.py -v
"""
import sys
import os
import warnings

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, HeatExchanger


@pytest.fixture(autouse=True)
def _quiet():
    warnings.filterwarnings('ignore')


def _pinch_tol(PPT):
    return max(1e-2, 1e-3 * PPT)


class TestEvaporatorHhiHhoThirdUnknown:
    """{'H_hi', 'H_ho', 'mh'} and {..., 'mc'}: only viable anchor is
    H_co -> H_hi, which used to demand Hot_Out.T up front."""

    def test_mh_unknown(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5)                      # H_hi, mh unknown
        M.add_point('water', 'ho', P=2e5)                      # H_ho unknown
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=10.0)
        M.add_point('water', 'co', P=3e5, T=340.0, Mass_flowrate=10.0)
        hx = HeatExchanger(M, 'E1', PPT=8.0, HEX_type='Evaporator',
                           HeatAdded=True,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
        assert hx.Hot_In.T > hx.Hot_Out.T
        assert hx.Hot_Mass_flowrate > 0

    def test_mc_unknown(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, Mass_flowrate=2.0)   # H_hi unknown
        M.add_point('water', 'ho', P=2e5, Mass_flowrate=2.0)   # H_ho unknown
        M.add_point('water', 'ci', P=3e5, T=300.0)
        M.add_point('water', 'co', P=3e5, T=340.0)             # mc unknown
        hx = HeatExchanger(M, 'E2', PPT=8.0, HEX_type='Evaporator',
                           HeatAdded=True,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
        assert hx.Hot_In.T > hx.Hot_Out.T
        assert hx.Cold_Mass_flowrate > 0


class TestCondenserHciHcoThirdUnknown:
    """{'H_ci', 'H_co', 'mh'} and {..., 'mc'}: only viable anchor is
    H_ho -> H_ci, which used to demand Cold_Out.T up front."""

    def test_mh_unknown(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, T=420.0)             # mh unknown
        M.add_point('water', 'ho', P=2e5, T=340.0)
        M.add_point('water', 'ci', P=3e5, Mass_flowrate=10.0)  # H_ci unknown
        M.add_point('water', 'co', P=3e5, Mass_flowrate=10.0)  # H_co unknown
        hx = HeatExchanger(M, 'C1', PPT=8.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
        assert hx.Cold_In.T < hx.Cold_Out.T
        assert hx.Hot_Mass_flowrate > 0

    def test_mc_unknown(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5, T=420.0, Mass_flowrate=2.0)
        M.add_point('water', 'ho', P=2e5, T=340.0, Mass_flowrate=2.0)
        M.add_point('water', 'ci', P=3e5)                      # H_ci unknown
        M.add_point('water', 'co', P=3e5)                      # H_co, mc unknown
        hx = HeatExchanger(M, 'C2', PPT=8.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
        assert hx.Cold_In.T < hx.Cold_Out.T
        assert hx.Cold_Mass_flowrate > 0


class TestDoublePipeUnaffectedAllAlong:
    """double_pipe always had the unaffected alternate anchor available
    for these same combinations -- confirm it still works post-fix."""

    def test_double_pipe_hi_ho_mh(self):
        M = ThermodynamicModel()
        M.add_point('water', 'hi', P=2e5)
        M.add_point('water', 'ho', P=2e5)
        M.add_point('water', 'ci', P=3e5, T=300.0, Mass_flowrate=10.0)
        M.add_point('water', 'co', P=3e5, T=340.0, Mass_flowrate=10.0)
        hx = HeatExchanger(M, 'D1', PPT=8.0, HEX_type='double_pipe',
                           HeatAdded=None,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.min_dT == pytest.approx(hx.PPT, abs=_pinch_tol(hx.PPT))
