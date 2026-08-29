"""
test_solve.py
--------------
Model.Solve() had no test coverage at all before this file: every other
test and example builds components in dependency order with
Calculate=True, bypassing the automatic solver entirely.

Run with:   pytest tests/test_solve.py -v
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import (
    ThermodynamicModel,
    Pump, Splitter, Mixer,
    HeatExchanger,
)


def fresh_model():
    m = ThermodynamicModel()
    m.set_dead_state(T0=298.15, P0=101325)
    return m


def _build_points(M):
    M.add_point('water', '1', P=5e6, T=350, Mass_flowrate=2.0)
    M.add_point('water', 'a', P=5e6)
    M.add_point('water', 'b', P=5e6)
    M.add_point('water', 'a2', P=6e6)
    M.add_point('water', 'b2', P=6e6)
    M.add_point('water', 'm', P=6e6)
    M.add_point('water', 'h_out', P=6e6, T=310)
    M.add_point('water', 'c1', P=2e5, T=290, Mass_flowrate=None)
    M.add_point('water', 'c2', P=2e5)
    M.add_point('water', 'h2_out', P=6e6, Mass_flowrate=None)
    M.add_point('water', 'd1', P=1e5, T=280, Mass_flowrate=1.5)
    M.add_point('water', 'd2', P=1e5, T=295)


def _build_reference():
    M = fresh_model()
    _build_points(M)
    Splitter(M, 'Splitter', '1', ['a', 'b'], [0.4, 0.6], Calculate=True)
    Pump(M, 'PumpA', 'a', 'a2', n_isen=0.8, Calculate=True)
    Pump(M, 'PumpB', 'b', 'b2', n_isen=0.8, Calculate=True)
    Mixer(M, 'Mixer', ['a2', 'b2'], 'm', Calculate=True)
    HeatExchanger(M, 'HEX1', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                  Hot_In_state='m', Hot_Out_state='h_out',
                  Cold_In_state='c1', Cold_Out_state='c2', Calculate=True)
    HeatExchanger(M, 'HEX2', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                  Hot_In_state='h_out', Hot_Out_state='h2_out',
                  Cold_In_state='d1', Cold_Out_state='d2', Calculate=True)
    return M


def _build_scrambled():
    M = fresh_model()
    _build_points(M)
    # Deliberately bad order: harder HEX first, pumps before the splitter
    # that feeds them, mixer before its inputs exist.
    HeatExchanger(M, 'HEX2', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                  Hot_In_state='h_out', Hot_Out_state='h2_out',
                  Cold_In_state='d1', Cold_Out_state='d2')
    Pump(M, 'PumpB', 'b', 'b2', n_isen=0.8)
    Mixer(M, 'Mixer', ['a2', 'b2'], 'm')
    HeatExchanger(M, 'HEX1', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                  Hot_In_state='m', Hot_Out_state='h_out',
                  Cold_In_state='c1', Cold_Out_state='c2')
    Pump(M, 'PumpA', 'a', 'a2', n_isen=0.8)
    Splitter(M, 'Splitter', '1', ['a', 'b'], [0.4, 0.6])
    ok = M.Solve(max_iter=10)
    assert ok is True
    return M


# HEX solves here are brentq-based (several _build_profile calls each, and
# each of those is itself expensive) -- module-scoped so the whole file
# pays for the reference and scrambled solves exactly once, not once per
# test.
@pytest.fixture(scope='module')
def ref():
    return _build_reference()


@pytest.fixture(scope='module')
def scr():
    return _build_scrambled()


class TestStagedSolveOrder:
    """
    Splitter -> two Pumps -> Mixer -> HEX1 (1-unknown pinch, mc) ->
    HEX2 (2-unknown pinch, mh + H_ho).

    Built two ways: in correct dependency order with Calculate=True
    (reference), and in a deliberately scrambled order -- the harder HEX
    registered before the easier one, pumps before their splitter, the
    mixer before its inputs exist -- solved with Model.Solve(). The
    default staged order (pass-through+junctions to convergence, then
    HEX ascending by current unknown-count, each tier to convergence)
    should reach the identical state regardless of construction order.
    """

    def test_scrambled_order_matches_reference(self, ref, scr):
        for name in ('1', 'a', 'b', 'a2', 'b2', 'm', 'h_out',
                     'c1', 'c2', 'h2_out', 'd1', 'd2'):
            h_ref = ref.Point[name].H
            h_scr = scr.Point[name].H
            assert h_ref is not None and h_scr is not None
            assert abs(h_ref - h_scr) <= 1e-3 * max(abs(h_ref), 1.0), name

    def test_scrambled_order_all_components_solved(self, scr):
        for cid, comp in scr.Component.items():
            assert comp.Solution_Status is True, cid

    def test_hex_unknowns_resolved(self, ref, scr):
        assert scr.Component['HEX1'].Cold_Mass_flowrate == pytest.approx(
            ref.Component['HEX1'].Cold_Mass_flowrate, rel=1e-6)
        assert scr.Component['HEX2'].Hot_Mass_flowrate == pytest.approx(
            ref.Component['HEX2'].Hot_Mass_flowrate, rel=1e-6)


class TestCustomOrder:
    """custom_order must still force the caller's exact solve order."""

    def _build_unordered(self, M):
        M.add_point('water', '1', P=5e6, T=350, Mass_flowrate=2.0)
        M.add_point('water', 'a', P=5e6)
        M.add_point('water', 'b', P=5e6)
        M.add_point('water', 'a2', P=6e6)
        M.add_point('water', 'b2', P=6e6)
        M.add_point('water', 'm', P=6e6)
        Splitter(M, 'Splitter', '1', ['a', 'b'], [0.4, 0.6])
        Pump(M, 'PumpA', 'a', 'a2', n_isen=0.8)
        Pump(M, 'PumpB', 'b', 'b2', n_isen=0.8)
        Mixer(M, 'Mixer', ['a2', 'b2'], 'm')

    def test_custom_order_still_solves(self):
        M = fresh_model()
        self._build_unordered(M)
        ok = M.Solve(custom_order=['Splitter', 'PumpA', 'PumpB', 'Mixer'])
        assert ok is True
        assert all(c.Solution_Status for c in M.Component.values())

    def test_custom_order_rejects_unknown_id(self):
        M = fresh_model()
        self._build_unordered(M)
        with pytest.raises(ValueError):
            M.Solve(custom_order=['NotAComponent'])
