"""
test_model.py
-------------
Tests for the ThermodynamicModel class itself
(add_point, Solve, ModelSummary, save/load).

Run with:   pytest tests/test_model.py -v
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import (
    ThermodynamicModel,
    Turbine, Pump, Splitter, Mixer, HeatExchanger,
)


def fresh_model():
    m = ThermodynamicModel()
    m.set_dead_state(T0=298.15, P0=101325)
    return m


class TestModelBasics:

    def test_add_point(self):
        M = fresh_model()
        M.add_point('water', '1', P=1e6, T=500)
        assert '1' in M.Point
        assert M.Point['1'].fluid == 'water'
        assert M.Point['1'].P == 1e6

    def test_multiple_points(self):
        M = fresh_model()
        M.add_point('water', '1', P=1e6, T=500)
        M.add_point('water', '2', P=2e6, T=600)
        assert len(M.Point) == 2

    def test_point_overwrite(self):
        """Adding same name overwrites the old point."""
        M = fresh_model()
        M.add_point('water', '1', P=1e6, T=500)
        old_H = M.Point['1'].H
        M.add_point('water', '1', P=1e6, T=600)
        assert M.Point['1'].H != old_H

    def test_str_no_crash(self):
        """str(Model) should not crash even when components aren't solved."""
        M = fresh_model()
        M.add_point('water', '1', P=1e6, T=500, Mass_flowrate=1)
        text = str(M)
        assert isinstance(text, str)


class TestModelSummary:

    def _build_simple_cycle(self):
        """Build a minimal Rankine cycle for testing."""
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
        M.add_point('water', '2', P=0.008e6)
        M.add_point('water', '3', P=0.008e6, Q=0)
        M.add_point('water', '4', P=8e6)

        Turbine(M, 'Turb', '1', '2', n_isen=0.85, Calculate=True)
        HeatExchanger(M, 'Cond', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=False,
                      Hot_In_state='2', Hot_Out_state='3',
                      Cold_In_state=None, Cold_Out_state=None,
                      Calculate=True)
        Pump(M, 'Pump', '3', '4', n_isen=1.0, Calculate=True)
        HeatExchanger(M, 'Boiler', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=True,
                      Hot_In_state=None, Hot_Out_state=None,
                      Cold_In_state='4', Cold_Out_state='1',
                      Calculate=True)
        return M

    def test_summary_runs(self):
        M = self._build_simple_cycle()
        M.ModelSummary()    # should not crash
        assert M.Net_power > 0
        assert M.Efficiency > 0
        assert M.Efficiency < 100

    def test_summary_twice(self):
        """Calling ModelSummary twice should not double values."""
        M = self._build_simple_cycle()
        M.ModelSummary()
        eff1 = M.Efficiency
        M.ModelSummary()
        eff2 = M.Efficiency
        assert abs(eff1 - eff2) < 0.01

    def test_point_print(self):
        M = self._build_simple_cycle()
        df = M.Point_print()
        assert len(df) == 4     # 4 state points


class TestSaveLoad:

    def _build_simple_cycle(self):
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
        M.add_point('water', '2', P=0.008e6)
        M.add_point('water', '3', P=0.008e6, Q=0)
        M.add_point('water', '4', P=8e6)

        Turbine(M, 'Turb', '1', '2', n_isen=0.85, Calculate=True)
        Pump(M, 'Pump', '3', '4', n_isen=1.0, Calculate=True)
        HeatExchanger(M, 'Cond', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=False,
                      Hot_In_state='2', Hot_Out_state='3',
                      Cold_In_state=None, Cold_Out_state=None,
                      Calculate=True)
        HeatExchanger(M, 'Boiler', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=True,
                      Hot_In_state=None, Hot_Out_state=None,
                      Cold_In_state='4', Cold_Out_state='1',
                      Calculate=True)
        return M

    def test_save_creates_file(self, tmp_path):
        """save_model should create a JSON file."""
        M = self._build_simple_cycle()
        filepath = os.path.join(str(tmp_path), 'test_save.json')
        M.save_model(filepath)
        assert os.path.isfile(filepath)

    def test_load_restores_points(self, tmp_path):
        """load_model should restore all state points."""
        M = self._build_simple_cycle()
        filepath = os.path.join(str(tmp_path), 'test_save.json')
        M.save_model(filepath)

        M2 = ThermodynamicModel()
        M2.load_model(filepath)

        assert '1' in M2.Point
        assert '2' in M2.Point
        assert '3' in M2.Point
        assert '4' in M2.Point

    def test_load_preserves_values(self, tmp_path):
        """Loaded point values should match original."""
        M = self._build_simple_cycle()
        filepath = os.path.join(str(tmp_path), 'test_save.json')
        M.save_model(filepath)

        M2 = ThermodynamicModel()
        M2.load_model(filepath)

        assert abs(M2.Point['1'].T - M.Point['1'].T) < 0.1
        assert abs(M2.Point['1'].P - M.Point['1'].P) < 1.0