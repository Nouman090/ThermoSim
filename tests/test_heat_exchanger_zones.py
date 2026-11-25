"""
Zone-by-zone UA for exchangers whose streams change phase.

A single LMTD from the four end temperatures assumes constant heat
capacities along the unit.  An evaporator, condenser or boiler is split
into phase zones instead (boiling, superheating, desuperheating,
condensing, ...), each with its own UA, and the total UA is their sum.
"""
import math
import warnings

import pytest
from CoolProp.CoolProp import PropsSI

from ThermoSim import ThermodynamicModel, HeatExchanger

F = 'n-Pentane'


def _model():
    M = ThermodynamicModel()
    M.set_dead_state()
    return M


def _evaporator(div_N=200):
    """R601 boiling at 5 C from x = 0.25 to 5 K superheat; water 14 -> 7 C."""
    M = _model()
    P = PropsSI('P', 'T', 278.15, 'Q', 1, F)
    M.add_point(F, 'e', P=P, Q=0.25)
    M.add_point(F, 'f', P=P, T=283.15)
    M.add_point('water', 'k1', T=287.15, P=3e5, Mass_flowrate=50.0)
    M.add_point('water', 'k2', T=280.15, P=3e5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        HeatExchanger(M, 'Evap', PPT=2, HEX_type='double_pipe', HeatAdded=None,
                      Hot_In_state='k1', Hot_Out_state='k2',
                      Cold_In_state='e', Cold_Out_state='f',
                      div_N=div_N, Calculate=True)
    return M, M.Component['Evap']


def _condenser():
    """R601 at 60 C, condensing at 40 C to saturated liquid; water 30 -> 35 C."""
    M = _model()
    P = PropsSI('P', 'T', 313.15, 'Q', 0, F)
    M.add_point(F, 'm', P=P, T=333.15, Mass_flowrate=40.0)
    M.add_point(F, 'd', P=P, Q=0)
    M.add_point('water', 'w1', T=303.15, P=3e5)
    M.add_point('water', 'w2', T=308.15, P=3e5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        HeatExchanger(M, 'Cond', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                      Hot_In_state='m', Hot_Out_state='d',
                      Cold_In_state='w1', Cold_Out_state='w2', Calculate=True)
    return M, M.Component['Cond']


def test_evaporator_has_boiling_and_superheating_zones():
    _, E = _evaporator()
    assert [(z['hot'], z['cold']) for z in E.zones] == \
        [('liquid', 'two-phase'), ('liquid', 'vapour')]
    assert list(E.UA_zones) == [z['label'] for z in E.zones]
    assert sum(z['Q'] for z in E.zones) == pytest.approx(E.Q, rel=1e-9)
    assert E.UA == pytest.approx(sum(E.UA_zones.values()), rel=1e-12)
    assert E.LMTD == pytest.approx(E.Q / E.UA, rel=1e-12)
    for z in E.zones:
        assert z['LMTD'] == pytest.approx(z['Q'] / z['UA'], rel=1e-12)


def test_boiling_zone_matches_hand_lmtd():
    """In the boiling zone the refrigerant is isothermal, so the zone UA is
    Q / LMTD of the water temperatures at the zone ends."""
    M, E = _evaporator()
    boil = E.zones[0]
    T_sat = PropsSI('T', 'P', M.Point['e'].P, 'Q', 1, F)
    d_hot_end = boil['T_hot_in'] - T_sat
    d_cold_end = boil['T_hot_out'] - T_sat
    lmtd = (d_hot_end - d_cold_end) / math.log(d_hot_end / d_cold_end)
    assert boil['UA'] == pytest.approx(boil['Q'] / lmtd, rel=2e-3)


def test_end_lmtd_is_wrong_for_phase_change_and_kept_separately():
    _, E = _evaporator()
    _, C = _condenser()
    # the single end-temperature LMTD overstates the evaporator UA and
    # understates the condenser UA here, by well over 10 %
    assert E.UA_ends > 1.3 * E.UA
    assert C.UA_ends < 0.9 * C.UA
    assert E.LMTD_ends == pytest.approx(E.Q / E.UA_ends, rel=1e-9)


def test_condenser_zones_desuperheating_then_condensing():
    _, C = _condenser()
    # listed from the hot-outlet / cold-inlet end
    assert [(z['hot'], z['cold']) for z in C.zones] == \
        [('two-phase', 'liquid'), ('vapour', 'liquid')]
    assert C.zones[0]['T_cold_in'] == pytest.approx(303.15, abs=1e-6)
    assert C.zones[-1]['T_hot_in'] == pytest.approx(333.15, abs=1e-6)
    assert C.zones[0]['T_hot_in'] == pytest.approx(C.zones[1]['T_hot_out'])


def test_single_phase_is_one_zone_close_to_end_lmtd():
    M = _model()
    M.add_point('water', 'a1', T=353.15, P=3e5, Mass_flowrate=10.0)
    M.add_point('water', 'a2', T=313.15, P=3e5)
    M.add_point('water', 'b1', T=293.15, P=3e5)
    M.add_point('water', 'b2', T=333.15, P=3e5)
    HeatExchanger(M, 'WW', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                  Hot_In_state='a1', Hot_Out_state='a2',
                  Cold_In_state='b1', Cold_Out_state='b2', Calculate=True)
    W = M.Component['WW']
    assert len(W.zones) == 1 and W.zones[0]['label'] == 'hot liquid, cold liquid'
    # only the variation of cp separates the two
    assert W.UA == pytest.approx(W.UA_ends, rel=2e-3)


def test_zone_ua_does_not_depend_on_div_N():
    _, coarse = _evaporator(div_N=20)
    _, fine = _evaporator(div_N=1000)
    assert coarse.UA == pytest.approx(fine.UA, rel=2e-3)
    for a, b in zip(coarse.zones, fine.zones):
        assert a['Q'] == pytest.approx(b['Q'], rel=1e-9)   # exact boundaries
        assert a['UA'] == pytest.approx(b['UA'], rel=2e-3)


def test_str_lists_the_zones():
    _, E = _evaporator()
    text = str(E)
    assert 'hot liquid, cold two-phase' in text
    assert 'hot liquid, cold vapour' in text


def test_oil_heated_boiler_has_three_zones():
    """Thermal oil (no two-phase envelope) heating R601 from subcooled
    liquid through boiling to superheat: preheating, boiling and
    superheating zones, with the oil labelled single-phase."""
    M = _model()
    P = PropsSI('P', 'T', 423.15, 'Q', 1, F)
    M.add_point('INCOMP::T66', 'h1', T=573.15, P=5e5, Mass_flowrate=20.0)
    M.add_point('INCOMP::T66', 'h2', P=5e5)
    M.add_point(F, 'a', P=P, T=313.15, Mass_flowrate=10.0)
    M.add_point(F, 'b', P=P, T=433.15)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        HeatExchanger(M, 'Gen', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                      Hot_In_state='h1', Hot_Out_state='h2',
                      Cold_In_state='a', Cold_Out_state='b', Calculate=True)
    G = M.Component['Gen']
    assert [(z['hot'], z['cold']) for z in G.zones] == [
        ('single-phase', 'liquid'), ('single-phase', 'two-phase'),
        ('single-phase', 'vapour')]
    assert G.UA == pytest.approx(sum(G.UA_zones.values()), rel=1e-12)
