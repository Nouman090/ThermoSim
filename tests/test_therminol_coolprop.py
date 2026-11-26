"""Therminol 66 comes from CoolProp's incompressible backend.

Up to v3.5.0 ThermoSim carried its own polynomial fit (T66.json) with an
enthalpy reference about 35 kJ/kg away from CoolProp's, so mixing ThermoSim
enthalpies with CoolProp temperatures put the oil ~17 K off.  Every
Therminol 66 name now resolves to 'INCOMP::T66' and is computed by CoolProp.
"""
import os
import warnings

import CoolProp.CoolProp as CP
import pytest

import ThermoSim
from ThermoSim import ThermodynamicModel, HeatExchanger
from ThermoSim.state import Prop, canonical_fluid, is_incompressible

T66 = 'INCOMP::T66'


@pytest.mark.parametrize('name', ['T66', 't66', 'Therminol66',
                                  'therminol-66', 'INCOMP::T66',
                                  'incomp::t66', ' T66 '])
def test_aliases_resolve_to_coolprop(name):
    assert canonical_fluid(name) == T66
    assert is_incompressible(name)
    p = Prop(name, 'x', T=450.0, P=5e5)
    assert p.fluid == T66


def test_other_fluids_unchanged():
    assert canonical_fluid('water') == 'water'
    assert canonical_fluid('R245fa') == 'R245fa'
    assert not is_incompressible('water')
    assert is_incompressible('INCOMP::Water')


@pytest.mark.parametrize('T', [300.0, 400.0, 500.0, 600.0])
def test_properties_equal_coolprop(T):
    p = Prop('T66', 'x', T=T, P=5e5)
    for attr, key in [('H', 'H'), ('S', 'S'), ('D', 'D'), ('Cp', 'C')]:
        assert getattr(p, attr) == pytest.approx(
            CP.PropsSI(key, 'T', T, 'P', 5e5, T66), rel=1e-12)


def test_no_quality_and_liquid_label():
    p = Prop('T66', 'x', T=450.0, P=5e5)
    assert p.Q is None
    assert p.phase == 'Incompressible liquid'
    assert p.Solution_Status is True


def test_enthalpy_input_round_trips():
    p1 = Prop('T66', 'a', T=523.15, P=5e5)
    p2 = Prop('T66', 'b', H=p1.H, P=5e5)
    assert p2.T == pytest.approx(523.15, abs=1e-6)
    assert p2.Q is None


def test_entropy_input_works():
    """Pumps set the isentropic outlet from (P, S)."""
    p1 = Prop('T66', 'a', T=423.15, P=2e5)
    p2 = Prop('T66', 'b', S=p1.S, P=10e5)
    assert p2.T == pytest.approx(423.15, abs=0.5)


def test_exergy_defined():
    M = ThermodynamicModel()
    p = M.add_point('T66', 'x', T=573.15, P=5e5, Mass_flowrate=10.0)
    assert p.ex is not None and p.ex > 0
    assert p.Ex == pytest.approx(10.0 * p.ex)


def test_update_keeps_alias_resolved():
    p = Prop('T66', 'x', T=450.0, P=5e5)
    p.update(T=500.0, P=5e5)
    assert p.fluid == T66
    assert p.H == pytest.approx(CP.PropsSI('H', 'T', 500, 'P', 5e5, T66))


def test_polynomial_data_file_removed():
    pkg = os.path.dirname(ThermoSim.__file__)
    assert not os.path.exists(os.path.join(pkg, 'T66.json'))


def test_heat_exchanger_with_oil_short_name():
    """An oil-heated boiler written with the short name solves, and the oil
    outlet temperature agrees with CoolProp at the computed enthalpy (the
    old polynomial was ~17 K off from CoolProp here)."""
    F = 'n-Pentane'
    M = ThermodynamicModel()
    M.set_dead_state()
    P = CP.PropsSI('P', 'T', 423.15, 'Q', 1, F)
    M.add_point('T66', 'h1', T=573.15, P=5e5, Mass_flowrate=20.0)
    M.add_point('T66', 'h2', P=5e5)
    M.add_point(F, 'a', P=P, T=313.15, Mass_flowrate=10.0)
    M.add_point(F, 'b', P=P, T=433.15)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        HeatExchanger(M, 'Gen', PPT=5, HEX_type='double_pipe', HeatAdded=None,
                      Hot_In_state='h1', Hot_Out_state='h2',
                      Cold_In_state='a', Cold_Out_state='b', Calculate=True)
    h1, h2 = M.Point['h1'], M.Point['h2']
    assert h2.fluid == T66
    assert h2.T == pytest.approx(
        CP.PropsSI('T', 'H', h2.H, 'P', h2.P, T66), abs=1e-6)
    a, b = M.Point['a'], M.Point['b']
    assert 20.0 * (h1.H - h2.H) == pytest.approx(
        10.0 * (b.H - a.H), rel=1e-9)
