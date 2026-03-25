"""
test_components.py
------------------
Tests for individual components.

Run with:   pytest tests/test_components.py -v
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Thermosim import (
    ThermodynamicModel,
    Turbine, Pump, Compressor,
    Mixer, Splitter, Expansion_valve,
    HeatExchanger,
    Source, Sink,
)


# ================================================================== #
#  Helper: create a fresh model with dead state
# ================================================================== #
def fresh_model():
    m = ThermodynamicModel()
    m.set_dead_state(T0=298.15, P0=101325)
    return m


# ================================================================== #
#  TURBINE TESTS
# ================================================================== #
class TestTurbine:

    def test_basic_expansion(self):
        """Turbine expands steam and produces positive work."""
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
        M.add_point('water', '2', P=0.7e6)

        Turbine(M, 'T1', '1', '2', n_isen=0.85, Calculate=True)

        assert M.Component['T1'].Solution_Status is True
        assert M.Component['T1'].work > 0       # turbine produces work
        assert M.Point['2'].T < M.Point['1'].T  # temperature drops

    def test_isentropic_turbine(self):
        """η=1 means outlet entropy equals inlet entropy."""
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
        M.add_point('water', '2', P=0.7e6)

        Turbine(M, 'T1', '1', '2', n_isen=1.0, Calculate=True)

        assert abs(M.Point['2'].S - M.Point['1'].S) < 1.0

    def test_mass_flow_propagation(self):
        """Mass flow should propagate from inlet to outlet."""
        M = fresh_model()
        M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=5.0)
        M.add_point('water', '2', P=0.7e6)

        Turbine(M, 'T1', '1', '2', n_isen=0.85, Calculate=True)

        assert M.Point['2'].Mass_flowrate == 5.0

    def test_work_scales_with_mass_flow(self):
        """Doubling mass flow should double work output."""
        M1 = fresh_model()
        M1.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
        M1.add_point('water', '2', P=0.7e6)
        Turbine(M1, 'T1', '1', '2', n_isen=0.85, Calculate=True)

        M2 = fresh_model()
        M2.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=2)
        M2.add_point('water', '2', P=0.7e6)
        Turbine(M2, 'T1', '1', '2', n_isen=0.85, Calculate=True)

        ratio = M2.Component['T1'].work / M1.Component['T1'].work
        assert abs(ratio - 2.0) < 0.01


# ================================================================== #
#  PUMP TESTS
# ================================================================== #
class TestPump:

    def test_basic_compression(self):
        """Pump compresses liquid and consumes work."""
        M = fresh_model()
        M.add_point('water', 'a', P=0.008e6, Q=0, Mass_flowrate=1)
        M.add_point('water', 'b', P=8e6)

        Pump(M, 'P1', 'a', 'b', n_isen=1.0, Calculate=True)

        assert M.Component['P1'].Solution_Status is True
        assert M.Component['P1'].work > 0       # pump consumes work
        assert M.Point['b'].P > M.Point['a'].P

    def test_isentropic_pump(self):
        """η=1 pump is isentropic."""
        M = fresh_model()
        M.add_point('water', 'a', P=0.008e6, Q=0, Mass_flowrate=1)
        M.add_point('water', 'b', P=8e6)

        Pump(M, 'P1', 'a', 'b', n_isen=1.0, Calculate=True)

        assert abs(M.Point['b'].S - M.Point['a'].S) < 1.0

    def test_lower_efficiency_more_work(self):
        """Lower efficiency → more work consumed."""
        M1 = fresh_model()
        M1.add_point('water', 'a', P=0.008e6, Q=0, Mass_flowrate=1)
        M1.add_point('water', 'b', P=8e6)
        Pump(M1, 'P1', 'a', 'b', n_isen=1.0, Calculate=True)

        M2 = fresh_model()
        M2.add_point('water', 'a', P=0.008e6, Q=0, Mass_flowrate=1)
        M2.add_point('water', 'b', P=8e6)
        Pump(M2, 'P1', 'a', 'b', n_isen=0.80, Calculate=True)

        assert M2.Component['P1'].work > M1.Component['P1'].work


# ================================================================== #
#  COMPRESSOR TESTS
# ================================================================== #
class TestCompressor:

    def test_basic_air_compression(self):
        """Compress air from 1 bar to 10 bar."""
        M = fresh_model()
        M.add_point('Air', 'c1', P=1e5, T=300, Mass_flowrate=1)
        M.add_point('Air', 'c2', P=10e5)

        Compressor(M, 'Comp1', 'c1', 'c2', n_isen=0.85, Calculate=True)

        comp = M.Component['Comp1']
        assert comp.Solution_Status is True
        assert comp.work > 0
        assert abs(comp.pressure_ratio - 10.0) < 0.01
        assert M.Point['c2'].T > M.Point['c1'].T

    def test_isentropic_compressor(self):
        """η=1 compressor is isentropic."""
        M = fresh_model()
        M.add_point('Air', 'c1', P=1e5, T=300, Mass_flowrate=1)
        M.add_point('Air', 'c2', P=10e5)

        Compressor(M, 'Comp1', 'c1', 'c2', n_isen=1.0, Calculate=True)

        assert abs(M.Point['c2'].S - M.Point['c1'].S) < 1.0

    def test_compressor_str(self):
        """__str__ should mention key info."""
        M = fresh_model()
        M.add_point('Air', 'c1', P=1e5, T=300, Mass_flowrate=1)
        M.add_point('Air', 'c2', P=10e5)
        Compressor(M, 'Comp1', 'c1', 'c2', n_isen=0.85, Calculate=True)

        text = str(M.Component['Comp1'])
        assert 'Compressor' in text
        assert 'Comp1' in text


# ================================================================== #
#  MIXER TESTS
# ================================================================== #
class TestMixer:

    def test_two_stream_mixing(self):
        """Mix two streams at the same pressure."""
        M = fresh_model()
        M.add_point('water', 's1', P=1e6, T=400, Mass_flowrate=2)
        M.add_point('water', 's2', P=1e6, T=350, Mass_flowrate=3)
        M.add_point('water', 'out', P=1e6)

        Mixer(M, 'Mix1', ['s1', 's2'], 'out', Calculate=True)

        assert M.Component['Mix1'].Solution_Status is True
        assert abs(M.Point['out'].Mass_flowrate - 5.0) < 0.01

    def test_energy_conservation(self):
        """Energy in = energy out for adiabatic mixer."""
        M = fresh_model()
        M.add_point('water', 's1', P=1e6, T=400, Mass_flowrate=2)
        M.add_point('water', 's2', P=1e6, T=350, Mass_flowrate=3)
        M.add_point('water', 'out', P=1e6)

        Mixer(M, 'Mix1', ['s1', 's2'], 'out', Calculate=True)

        E_in = (M.Point['s1'].H * M.Point['s1'].Mass_flowrate +
                M.Point['s2'].H * M.Point['s2'].Mass_flowrate)
        E_out = M.Point['out'].H * M.Point['out'].Mass_flowrate
        assert abs(E_in - E_out) < 1.0


# ================================================================== #
#  SPLITTER TESTS
# ================================================================== #
class TestSplitter:

    def test_split_50_50(self):
        """Split a stream 50/50."""
        M = fresh_model()
        M.add_point('water', 'in', P=1e6, T=500, Mass_flowrate=10)
        M.add_point('water', 'a', P=1e6)
        M.add_point('water', 'b', P=1e6)

        Splitter(M, 'Sp1', 'in', ['a', 'b'], [0.5, 0.5], Calculate=True)

        assert abs(M.Point['a'].Mass_flowrate - 5.0) < 0.01
        assert abs(M.Point['b'].Mass_flowrate - 5.0) < 0.01
        # enthalpy preserved
        assert abs(M.Point['a'].H - M.Point['in'].H) < 0.01

    def test_split_mass_conservation(self):
        """Total outlet mass = inlet mass."""
        M = fresh_model()
        M.add_point('water', 'in', P=1e6, T=500, Mass_flowrate=10)
        M.add_point('water', 'a', P=1e6)
        M.add_point('water', 'b', P=1e6)
        M.add_point('water', 'c', P=1e6)

        Splitter(M, 'Sp1', 'in', ['a', 'b', 'c'],
                 [0.5, 0.3, 0.2], Calculate=True)

        total = sum(M.Point[n].Mass_flowrate for n in ['a', 'b', 'c'])
        assert abs(total - 10.0) < 0.01


# ================================================================== #
#  EXPANSION VALVE TESTS
# ================================================================== #
class TestExpansionValve:

    def test_isenthalpic(self):
        """Enthalpy is conserved across expansion valve."""
        M = fresh_model()
        M.add_point('R134a', 'v1', P=1e6, T=310, Mass_flowrate=1)
        M.add_point('R134a', 'v2', P=2e5)

        Expansion_valve(M, 'EV1', 'v1', 'v2', Calculate=True)

        assert abs(M.Point['v2'].H - M.Point['v1'].H) < 1.0


# ================================================================== #
#  HEAT EXCHANGER TESTS  (SimpleHEX)
# ================================================================== #
class TestSimpleHEX:

    def test_cold_side_heating(self):
        """SimpleHEX: known Q heats the cold side."""
        M = fresh_model()
        M.add_point('water', 'ci', P=1e6, T=350, Mass_flowrate=1)
        M.add_point('water', 'co', P=1e6)

        HeatExchanger(M, 'Heater', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=True,
                      Hot_In_state=None, Hot_Out_state=None,
                      Cold_In_state='ci', Cold_Out_state='co',
                      Q=500e3, Calculate=True)

        assert M.Component['Heater'].Solution_Status is True
        assert M.Point['co'].T > M.Point['ci'].T

    def test_hot_side_cooling(self):
        """SimpleHEX: known Q cools the hot side."""
        M = fresh_model()
        M.add_point('water', 'hi', P=1e6, T=500, Mass_flowrate=1)
        M.add_point('water', 'ho', P=1e6)

        HeatExchanger(M, 'Cooler', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=False,
                      Hot_In_state='hi', Hot_Out_state='ho',
                      Cold_In_state=None, Cold_Out_state=None,
                      Q=200e3, Calculate=True)

        assert M.Point['ho'].T < M.Point['hi'].T

    def test_Q_calculated(self):
        """SimpleHEX: Q computed when both sides known."""
        M = fresh_model()
        M.add_point('water', 'hi', P=1e6, T=500, Mass_flowrate=1)
        M.add_point('water', 'ho', P=1e6, T=400)

        HeatExchanger(M, 'Calc', PPT=5, HEX_type='SimpleHEX',
                      HeatAdded=False,
                      Hot_In_state='hi', Hot_Out_state='ho',
                      Cold_In_state=None, Cold_Out_state=None,
                      Calculate=True)

        assert M.Component['Calc'].Q is not None
        assert M.Component['Calc'].Q > 0


# ================================================================== #
#  SOURCE / SINK TESTS
# ================================================================== #
class TestSourceSink:

    def test_source(self):
        M = fresh_model()
        M.add_point('water', 'out', P=1e6, T=500, Mass_flowrate=1)
        Source(M, 'Src', 'out', Calculate=True)
        assert M.Component['Src'].Solution_Status is True
        assert M.Component['Src'].energy_supply > 0

    def test_sink(self):
        M = fresh_model()
        M.add_point('water', 'in', P=1e6, T=500, Mass_flowrate=1)
        Sink(M, 'Snk', 'in', Calculate=True)
        assert M.Component['Snk'].Solution_Status is True