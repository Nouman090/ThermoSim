"""
test_state.py
-------------
Tests for the Prop (state point) class.

Run with:   pytest tests/test_state.py -v
"""

import pytest
import sys
import os

# Make sure the project root is in the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel
from ThermoSim.state import Prop
from ThermoSim import config


class TestPropCreation:
    """Test that Prop objects are created correctly."""

    def setup_method(self):
        """Run before each test."""
        config.set_dead_state(T0=298.15, P0=101325)

    def test_water_PT(self):
        """Create a water state from P and T."""
        pt = Prop('water', 'test_PT', P=101325, T=373.15)
        assert pt.Solution_Status is True
        assert pt.P == 101325
        assert pt.T == 373.15
        assert pt.H is not None
        assert pt.S is not None
        assert pt.fluid == 'water'

    def test_water_PH(self):
        """Create a water state from P and H."""
        # first get H at known conditions
        ref = Prop('water', 'ref', P=1e6, T=500)
        pt = Prop('water', 'test_PH', P=1e6, H=ref.H)
        assert pt.Solution_Status is True
        assert abs(pt.T - 500) < 0.1

    def test_water_PS(self):
        """Create a water state from P and S."""
        ref = Prop('water', 'ref', P=1e6, T=500)
        pt = Prop('water', 'test_PS', P=1e6, S=ref.S)
        assert pt.Solution_Status is True
        assert abs(pt.H - ref.H) < 1.0

    def test_water_PQ_saturated_liquid(self):
        """Saturated liquid at given pressure."""
        pt = Prop('water', 'sat_liq', P=101325, Q=0)
        assert pt.Solution_Status is True
        assert abs(pt.T - 373.15) < 0.5   # ~100°C at 1 atm

    def test_water_PQ_saturated_vapor(self):
        """Saturated vapor at given pressure."""
        pt = Prop('water', 'sat_vap', P=101325, Q=1)
        assert pt.Solution_Status is True
        assert abs(pt.T - 373.15) < 0.5

    def test_mass_flowrate(self):
        """Mass flow rate is stored correctly."""
        pt = Prop('water', 'test_mf', P=1e6, T=500, Mass_flowrate=2.5)
        assert pt.Mass_flowrate == 2.5

    def test_exergy_calculated(self):
        """Exergy is calculated when mass flow is given."""
        pt = Prop('water', 'test_ex', P=1e6, T=500, Mass_flowrate=1.0)
        assert pt.ex is not None
        assert pt.Ex is not None
        assert pt.Ex == pt.ex * pt.Mass_flowrate

    def test_exergy_none_without_massflow(self):
        """Ex should be None when mass flow is not given."""
        pt = Prop('water', 'test_ex2', P=1e6, T=500)
        assert pt.ex is not None     # specific exergy exists
        assert pt.Ex is None         # total exergy needs mass flow

    def test_invalid_property_raises(self):
        """Passing an invalid property key should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid property"):
            Prop('water', 'bad', P=1e6, X=500)

    def test_too_many_properties_raises(self):
        """More than 2 properties should raise ValueError."""
        with pytest.raises(ValueError, match="exactly 2"):
            Prop('water', 'toomany', P=1e6, T=500, H=3000e3)

    def test_different_fluids(self):
        """Test with R134a (a refrigerant)."""
        pt = Prop('R134a', 'r134a_test', P=500e3, T=300)
        assert pt.Solution_Status is True
        assert pt.H is not None

    def test_str_method(self):
        """Test that __str__ doesn't crash."""
        pt = Prop('water', 'str_test', P=1e6, T=500, Mass_flowrate=1.0)
        text = str(pt)
        assert 'water' in text
        assert 'str_test' in text

    def test_repr_method(self):
        """Test that __repr__ doesn't crash."""
        pt = Prop('water', 'repr_test', P=1e6, T=500)
        text = repr(pt)
        assert 'repr_test' in text