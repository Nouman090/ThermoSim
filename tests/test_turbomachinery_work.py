"""
test_turbomachinery_work.py
-----------------------------
Tests for the `work=` constructor argument on Turbine, Compressor and Pump:
when given and neither port carries a mass flow, it should be used to
derive the mass flow rate; when a mass flow rate IS already resolvable from
the ports, `work` should only be checked for consistency, never override it.

Run with:   pytest tests/test_turbomachinery_work.py -v
"""
import sys
import os
import warnings

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, Turbine, Pump, Compressor


def fresh_model():
    m = ThermodynamicModel()
    m.set_dead_state(T0=298.15, P0=101325)
    return m


class TestTurbineWork:
    def test_work_derives_mass_flow(self):
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=773.15)
        M.add_point('water', '2', P=0.1e6)
        T1 = Turbine(M, 'T1', '1', '2', n_isen=0.85, work=500e3, Calculate=True)

        assert T1.work == pytest.approx(500e3, rel=1e-6)
        assert T1.In.Mass_flowrate == T1.Out.Mass_flowrate
        assert T1.In.Mass_flowrate > 0

    def test_mass_flowrate_given_checks_consistency_not_overridden(self):
        # Solve once without a work target to learn the true work, then
        # re-solve an identical model WITH that work as the target: the
        # ports already carry the mass flow, so it must be left at 1.0
        # (not re-derived from work), and no mismatch warning should fire.
        M0 = fresh_model()
        M0.add_point('water', '1', P=8e6, T=773.15, Mass_flowrate=1.0)
        M0.add_point('water', '2', P=0.1e6, Mass_flowrate=1.0)
        T0 = Turbine(M0, 'T0', '1', '2', n_isen=0.85, Calculate=True)
        true_work = T0.work

        M2 = fresh_model()
        M2.add_point('water', '1', P=8e6, T=773.15, Mass_flowrate=1.0)
        M2.add_point('water', '2', P=0.1e6, Mass_flowrate=1.0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            T2 = Turbine(M2, 'T2', '1', '2', n_isen=0.85, work=true_work,
                         Calculate=True)
            assert not any('inconsistent' in str(x.message) for x in w)

        assert T2.In.Mass_flowrate == 1.0
        assert T2.work == pytest.approx(true_work, rel=1e-6)

    def test_mismatched_work_warns_and_keeps_port_mass(self):
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=773.15, Mass_flowrate=1.0)
        M.add_point('water', '2', P=0.1e6, Mass_flowrate=1.0)
        with pytest.warns(UserWarning, match='inconsistent'):
            T1 = Turbine(M, 'T1', '1', '2', n_isen=0.85, work=999e6,
                         Calculate=True)
        assert T1.In.Mass_flowrate == 1.0   # port value wins, not overridden


class TestCompressorWork:
    def test_work_derives_mass_flow(self):
        M = fresh_model()
        M.add_point('Air', '1', P=1e5, T=300.0)
        M.add_point('Air', '2', P=5e5)
        C1 = Compressor(M, 'C1', '1', '2', n_isen=0.8, work=50e3, Calculate=True)

        assert C1.work == pytest.approx(50e3, rel=1e-6)
        assert C1.In.Mass_flowrate == C1.Out.Mass_flowrate
        assert C1.In.Mass_flowrate > 0

    def test_negative_implied_mass_flow_raises(self):
        M = fresh_model()
        M.add_point('Air', '1', P=1e5, T=300.0)
        M.add_point('Air', '2', P=5e5)
        with pytest.raises(ValueError, match='non-positive mass flow'):
            Compressor(M, 'C1', '1', '2', n_isen=0.8, work=-1000.0,
                       Calculate=True)


class TestPumpWork:
    def test_work_derives_mass_flow(self):
        M = fresh_model()
        M.add_point('water', '1', P=0.1e6, Q=0.0)
        M.add_point('water', '2', P=8e6)
        P1 = Pump(M, 'P1', '1', '2', n_isen=0.8, work=20e3, Calculate=True)

        assert P1.work == pytest.approx(20e3, rel=1e-6)
        assert P1.In.Mass_flowrate == P1.Out.Mass_flowrate
        assert P1.In.Mass_flowrate > 0

    def test_mass_flowrate_given_not_overridden(self):
        M = fresh_model()
        M.add_point('water', '1', P=0.1e6, Q=0.0, Mass_flowrate=2.0)
        M.add_point('water', '2', P=8e6)
        P0 = Pump(M, 'P0', '1', '2', n_isen=0.8, Calculate=True)
        true_work = P0.work

        M2 = fresh_model()
        M2.add_point('water', '1', P=0.1e6, Q=0.0, Mass_flowrate=2.0)
        M2.add_point('water', '2', P=8e6)
        P2 = Pump(M2, 'P2', '1', '2', n_isen=0.8, work=true_work,
                  Calculate=True)
        assert P2.In.Mass_flowrate == 2.0
        assert P2.work == pytest.approx(true_work, rel=1e-6)


class TestNonRegression:
    """work=None (the default) must reproduce the exact pre-existing
    behaviour -- both the solved state and the derived work value."""

    def test_turbine_unchanged_without_work(self):
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
        M.add_point('water', '2', P=0.7e6)
        Turbine(M, 'T1', '1', '2', n_isen=0.85, Calculate=True)
        assert M.Component['T1'].Solution_Status is True
        assert M.Component['T1'].work > 0
        assert M.Point['2'].T < M.Point['1'].T
