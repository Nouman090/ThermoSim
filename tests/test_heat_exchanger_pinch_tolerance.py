"""
test_heat_exchanger_pinch_tolerance.py
--------------------------------------
Regression test for a false "PPT violation" from HeatExchanger.

_solve_2_unknown located the pinch with brentq(xtol=1e-4, rtol=1e-4) on an
enthalpy near 4.4e5 J/kg, so the search could stop some tens of J/kg from the
root, a few hundredths of a kelvin.  _check_PPT_ends allows only
max(1e-2, 1e-3*PPT) K, which is 0.01 K at a 5 K pinch, and rejected the result:
an R1233zd(E) regenerator failed at PPT = 4.9 to 5.5 K while its neighbours
passed, depending on where brentq happened to stop.

The fix tightens the search rather than widening the allowance, so these tests
check both directions: good solutions land on the target, and real violations
are still refused.

Run with:   pytest tests/test_heat_exchanger_pinch_tolerance.py -v
"""
import sys
import os
import warnings

import pytest
from CoolProp.CoolProp import PropsSI

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, HeatExchanger


FLUID = 'R1233zd(E)'
P_HOT = PropsSI('P', 'T', 313.15, 'Q', 0, FLUID)
P_COLD = PropsSI('P', 'T', 368.3, 'Q', 1, FLUID)

#: Far below the end-check allowance, and above what the old tolerances reached
#: on this case (0.014 K at PPT = 5 K).
PINCH_ACCURACY = 1e-4


@pytest.fixture(autouse=True)
def _quiet():
    warnings.filterwarnings('ignore')


def _exchanger(model, ppt):
    return HeatExchanger(model, 'HX1', PPT=ppt, HEX_type='double_pipe',
                         HeatAdded=True,
                         Hot_In_state='hi', Hot_Out_state='ho',
                         Cold_In_state='ci', Cold_Out_state='co')


def _regenerator(ppt):
    """Reported case: both outlets unknown, P and mass flow given."""
    M = ThermodynamicModel()
    M.add_point(FLUID, 'hi', P=P_HOT, T=322.75, Mass_flowrate=1.838)
    M.add_point(FLUID, 'ho', P=P_HOT, Mass_flowrate=1.838)
    M.add_point(FLUID, 'ci', P=P_COLD, T=313.52, Mass_flowrate=1.0)
    M.add_point(FLUID, 'co', P=P_COLD, Mass_flowrate=1.0)
    return M, _exchanger(M, ppt)


class TestReportedRegenerator:
    # 4.9 to 5.5 K all failed before the fix; 4 K and 6 K passed the check but
    # sat 0.012 K and 0.008 K from the target.
    @pytest.mark.parametrize('ppt', [1, 4, 4.9, 5, 5.03, 5.5, 6, 9])
    def test_solves_on_the_target_pinch(self, ppt):
        _, hx = _regenerator(ppt)
        hx.Cal()
        assert hx.min_dT == pytest.approx(ppt, abs=PINCH_ACCURACY)
        # The pinch of this exchanger sits at the hot-inlet end.
        assert hx.Hot_In.T - hx.Cold_Out.T == pytest.approx(ppt, abs=PINCH_ACCURACY)
        assert hx.Hot_Out.T - hx.Cold_In.T > ppt

    def test_energy_balance_closes(self):
        _, hx = _regenerator(5.0)
        hx.Cal()
        Q_hot = hx.Hot_Mass_flowrate * (hx.Hot_In.H - hx.Hot_Out.H)
        Q_cold = hx.Cold_Mass_flowrate * (hx.Cold_Out.H - hx.Cold_In.H)
        assert Q_hot == pytest.approx(Q_cold, rel=1e-9)
        assert Q_hot > 0


class TestOtherSearchesAreTight:
    """Every pinch search shares the tolerances, not only {H_ho, H_co}."""

    def _points(self, M, known):
        full = {'hi': 400.0, 'ho': 340.0, 'ci': 300.0, 'co': 360.0}
        for name, T in full.items():
            flow = 1.0 if ('mh' if name[0] == 'h' else 'mc') in known else None
            kwargs = {'T': T} if name in known else {}
            M.add_point('Water', name, P=10e5, Mass_flowrate=flow, **kwargs)

    @pytest.mark.parametrize('known', [
        ('hi', 'ci', 'mh', 'mc'),    # {H_ho, H_co}
        ('ho', 'ci', 'co', 'mh'),    # {H_hi, mc}
        ('ci', 'co', 'mh', 'mc'),    # {H_hi, H_ho}
        ('hi', 'ho', 'mh', 'mc'),    # {H_ci, H_co}
        ('hi', 'ho', 'ci', 'mh'),    # {H_co, mc}
        ('ci', 'co', 'hi', 'mc'),    # {H_ho, mh}
    ])
    def test_pinch_is_hit(self, known):
        M = ThermodynamicModel()
        self._points(M, known)
        hx = _exchanger(M, 7.0)
        hx.Cal()
        assert hx.min_dT == pytest.approx(7.0, abs=PINCH_ACCURACY)


class TestRealViolationsStillRejected:
    """The allowance itself did not widen."""

    def _fully_known(self, dT_hot_end, ppt=5.0):
        M = ThermodynamicModel()
        M.add_point('Water', 'hi', P=1e5, T=350.0, Mass_flowrate=1.0)
        M.add_point('Water', 'ho', P=1e5, Mass_flowrate=1.0)
        M.add_point('Water', 'ci', P=1e5, T=300.0, Mass_flowrate=1.0)
        M.add_point('Water', 'co', P=1e5, T=350.0 - dT_hot_end, Mass_flowrate=1.0)
        return _exchanger(M, ppt)

    @pytest.mark.parametrize('dT', [1.0, 4.9, 4.98])
    def test_end_violation_raises(self, dT):
        with pytest.raises(ValueError, match='PPT violation'):
            self._fully_known(dT).Cal()

    @pytest.mark.parametrize('dT', [4.995, 5.0, 5.1])
    def test_within_allowance_passes(self, dT):
        hx = self._fully_known(dT)
        hx.Cal()
        assert hx.Solution_Status is True


class TestInternalPinch:
    """min_dT still comes from the whole profile, not only the ends."""

    def test_boiling_stream_pinches_inside(self):
        P_boil = PropsSI('P', 'T', 370.0, 'Q', 0, 'R245fa')
        M = ThermodynamicModel()
        M.add_point('Water', 'hi', P=5e5, T=403.0, Mass_flowrate=1.0)
        M.add_point('Water', 'ho', P=5e5, Mass_flowrate=1.0)
        M.add_point('R245fa', 'ci', P=P_boil, T=310.0, Mass_flowrate=1.0)
        M.add_point('R245fa', 'co', P=P_boil, Mass_flowrate=1.0)
        hx = _exchanger(M, 8.0)
        hx.Cal()
        assert hx.min_dT == pytest.approx(8.0, abs=PINCH_ACCURACY)
        assert hx.Hot_In.T - hx.Cold_Out.T > 8.0 + 1.0
        assert hx.Hot_Out.T - hx.Cold_In.T > 8.0 + 1.0
