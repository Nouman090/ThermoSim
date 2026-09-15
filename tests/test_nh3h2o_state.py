"""
test_nh3h2o_state.py
---------------------
Tests for the NH3-H2O mixture state layer: Prop._calc_nh3h2o and
NH3H2OFluid (state.py, nh3h2o.py).

These mirror NH3H2O/validate_lib.py's own checks, but run them through
Prop/NH3H2OFluid rather than the bare library, so the SI unit-conversion
layer (Pa<->bar, J/kg<->kJ/kg, J/kg/K<->kJ/kg/K) is what's actually being
verified -- nh3h2o.py's own numerics are already validated standalone.

Run with:   pytest tests/test_nh3h2o_state.py -v
"""
import sys
import os
import warnings

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, NH3H2OFluid
from ThermoSim import config
import ThermoSim.nh3h2o as aw


@pytest.fixture(autouse=True)
def _dead_state():
    config.set_dead_state(T0=298.15, P0=101325)


class TestNH3H2OFluid:
    def test_valid_composition(self):
        f = NH3H2OFluid(0.4)
        assert f.z == pytest.approx(0.4)

    @pytest.mark.parametrize("z", [-0.1, 1.1])
    def test_invalid_composition_raises(self, z):
        with pytest.raises(ValueError):
            NH3H2OFluid(z)

    def test_equality_and_hash(self):
        assert NH3H2OFluid(0.4) == NH3H2OFluid(0.4)
        assert NH3H2OFluid(0.4) != NH3H2OFluid(0.5)
        assert hash(NH3H2OFluid(0.4)) == hash(NH3H2OFluid(0.4))

    def test_repr(self):
        assert 'z=0.4000' in repr(NH3H2OFluid(0.4))


class TestPropDispatch:
    """Prop must route NH3H2OFluid states through _calc_nh3h2o, not CoolProp."""

    def test_is_nh3h2o_dispatch(self):
        M = ThermodynamicModel()
        pt = M.add_point(NH3H2OFluid(0.4), 'pt', P=4e5, Q=0.0)
        assert pt.Solution_Status is True
        assert pt.z == pytest.approx(0.4)


class TestUnitConversionRoundTrip:
    """Every Prop-computed property must match nh3h2o.py's own answer after
    converting units by hand (Pa->bar, J/kg->kJ/kg, J/kg/K->kJ/kg/K)."""

    @pytest.mark.parametrize("z,P_bar,T", [
        (0.1, 4.0, 300.0), (0.4, 15.0, 400.0), (0.7, 10.0, 280.0),
        (0.9999, 8.0, 320.0),  # near-pure ammonia
        (0.0001, 8.0, 380.0),  # near-pure water
    ])
    def test_PT_matches_reference(self, z, P_bar, T):
        ref = aw.state(z=z, P=P_bar, T=T, warn=False)
        M = ThermodynamicModel()
        pt = M.add_point(NH3H2OFluid(z), 'pt', P=P_bar * 1e5, T=T)
        assert pt.T == pytest.approx(ref.T, abs=1e-6)
        assert pt.P == pytest.approx(ref.P * 1e5, rel=1e-9)
        assert pt.H == pytest.approx(ref.h * 1e3, rel=1e-9)
        assert pt.S == pytest.approx(ref.s * 1e3, rel=1e-9)
        assert pt.D == pytest.approx(1.0 / ref.v, rel=1e-9)

    def test_bubble_point_PQ(self):
        z, P_bar = 0.4, 4.0
        ref = aw.state(z=z, P=P_bar, Q=0.0, warn=False)
        M = ThermodynamicModel()
        pt = M.add_point(NH3H2OFluid(z), 'pt', P=P_bar * 1e5, Q=0.0)
        assert pt.T == pytest.approx(ref.T, abs=1e-6)
        assert pt.H == pytest.approx(ref.h * 1e3, rel=1e-9)
        assert pt.Q == pytest.approx(0.0, abs=1e-9)

    def test_two_phase_flash_populates_xy(self):
        z, P_bar, T = 0.4, 15.0, 400.0
        ref = aw.state(z=z, P=P_bar, T=T, warn=False)
        M = ThermodynamicModel()
        pt = M.add_point(NH3H2OFluid(z), 'pt', P=P_bar * 1e5, T=T)
        assert ref.phase == 'two-phase'   # sanity on the chosen fixture
        assert pt.x == pytest.approx(ref.x, abs=1e-6)
        assert pt.y == pytest.approx(ref.y, abs=1e-6)
        assert pt.Q == pytest.approx(ref.Q, abs=1e-6)
        assert 'Two-phase' in pt.phase

    def test_PH_round_trip(self):
        z = 0.4
        M = ThermodynamicModel()
        pt = M.add_point(NH3H2OFluid(z), 'pt', P=4e5, Q=0.0)
        pt2 = M.add_point(NH3H2OFluid(z), 'pt2', P=pt.P, H=pt.H)
        assert pt2.T == pytest.approx(pt.T, rel=1e-6)


class TestAmmoniaWaterErrorMapping:
    def test_invalid_pair_raises_value_error(self):
        M = ThermodynamicModel()
        with pytest.raises(ValueError):
            # T,H in the compressed liquid is a documented non-unique pair;
            # this particular combination has no solution at all, which is
            # the simplest reliable way to trigger AmmoniaWaterError.
            M.add_point(NH3H2OFluid(0.4), 'bad', T=230.0, H=-1e9)


class TestExergy:
    def test_zero_at_dead_state(self):
        M = ThermodynamicModel()
        M.set_dead_state(T0=298.15, P0=101325)
        pt = M.add_point(NH3H2OFluid(0.4), 'pt',
                         T=298.15, P=101325, Mass_flowrate=1.0)
        assert pt.ex == pytest.approx(0.0, abs=1e-2)
        assert pt.Ex == pytest.approx(0.0, abs=1e-2)

    def test_mass_flowrate_reassignment_recomputes_Ex(self):
        M = ThermodynamicModel()
        M.set_dead_state()
        pt = M.add_point(NH3H2OFluid(0.4), 'pt', P=4e5, Q=0.0, Mass_flowrate=1.0)
        assert pt.Ex is not None
        pt.Mass_flowrate = 2.0
        assert pt.Ex == pytest.approx(pt.ex * 2.0, rel=1e-9)


class TestNonRegression:
    """Adding the NH3H2O dispatch branch must not disturb CoolProp/Therminol
    fluids -- covered exhaustively by the rest of tests/, this is a quick
    sanity check local to this file."""

    def test_water_unaffected(self):
        M = ThermodynamicModel()
        pt = M.add_point('water', 'w', P=101325, T=373.15)
        assert pt.fluid == 'water'
        assert pt.Solution_Status is True
