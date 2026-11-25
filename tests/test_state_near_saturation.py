"""
test_state_near_saturation.py
-----------------------------
Regression tests for the near-saturation snap in ``Prop._calc_coolprop``.

CoolProp refuses a (T, P) flash that sits within ~1e-4 % of the saturation
pressure, so ThermoSim snaps (T, P) inputs within a small relative band of
T_sat onto the saturation line.  Up to v3.5.0 the snap always picked
saturated LIQUID (Q=0), even for a state slightly ABOVE T_sat: a suction
state 0.022 K superheated came back as liquid with h_l instead of h_v,
and a compressor fed from it produced an outlet 34 K below saturation.

The snap is now side-aware: T above T_sat -> saturated vapour (Q=1),
T below (or equal to, within round-off) T_sat -> saturated liquid (Q=0).

Run with:   pytest tests/test_state_near_saturation.py -v
"""

import os
import sys
import warnings

import pytest
import CoolProp.CoolProp as CP

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, config
from ThermoSim.state import Prop


PENTANE = 'n-Pentane'
T10 = 279.15                                   # K, the thesis evaporator
P10 = CP.PropsSI('P', 'T', T10, 'Q', 0, PENTANE)
H_LIQ = CP.PropsSI('H', 'P', P10, 'Q', 0, PENTANE)
H_VAP = CP.PropsSI('H', 'P', P10, 'Q', 1, PENTANE)


def _make(fluid, T, P, name='x'):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        pt = Prop(fluid, name, T=T, P=P)
    msgs = [str(x.message) for x in w if 'near saturation' in str(x.message)]
    return pt, msgs


class TestPureFluidNearSaturation:

    def setup_method(self):
        config.set_dead_state(T0=298.15, P0=101325)

    def test_slightly_superheated_snaps_to_vapour(self):
        """The reported case: 0.022 K above T_sat must be vapour, not liquid."""
        pt, msgs = _make(PENTANE, T10 + 0.022, P10)
        assert pt.Solution_Status is True
        assert pt.phase == 'Saturated vapour'
        assert pt.Q == 1.0
        # Saturated vapour enthalpy; direct CoolProp gives ~311441 J/kg,
        # so the snap is within cp*dT (~40 J/kg) of the true state.
        assert pt.H == pytest.approx(H_VAP, rel=1e-9)
        direct = CP.PropsSI('H', 'T', T10 + 0.022, 'P', P10, PENTANE)
        assert abs(pt.H - direct) < 100.0
        assert pt.H > H_LIQ + 3e5          # nowhere near the liquid value
        assert msgs and 'saturated vapour (Q=1)' in msgs[0]

    def test_slightly_subcooled_snaps_to_liquid(self):
        """Below T_sat the snap is unchanged: saturated liquid."""
        pt, msgs = _make(PENTANE, T10 - 0.022, P10)
        assert pt.phase == 'Saturated liquid'
        assert pt.Q == 0.0
        assert pt.H == pytest.approx(H_LIQ, rel=1e-9)
        assert msgs and 'saturated liquid (Q=0)' in msgs[0]

    def test_exactly_at_saturation_stays_liquid(self):
        """T == T_sat(P) (round-off only) keeps the historical Q=0 choice."""
        pt, msgs = _make(PENTANE, T10, P10)
        assert pt.phase == 'Saturated liquid'
        assert pt.H == pytest.approx(H_LIQ, rel=1e-9)
        assert msgs

    @pytest.mark.parametrize('dT', [-0.5, 0.5])
    def test_outside_band_is_not_snapped(self, dT):
        """Beyond the band, CoolProp resolves the state directly."""
        pt, msgs = _make(PENTANE, T10 + dT, P10)
        assert not msgs
        assert pt.Q is None
        assert pt.H == pytest.approx(
            CP.PropsSI('H', 'T', T10 + dT, 'P', P10, PENTANE), rel=1e-9)

    def test_water_both_sides(self):
        P = 101325.0
        Tsat = CP.PropsSI('T', 'P', P, 'Q', 0, 'water')
        above, _ = _make('water', Tsat + 0.02, P)
        below, _ = _make('water', Tsat - 0.02, P)
        assert above.phase == 'Saturated vapour'
        assert below.phase == 'Saturated liquid'
        assert above.H - below.H == pytest.approx(
            CP.PropsSI('H', 'P', P, 'Q', 1, 'water')
            - CP.PropsSI('H', 'P', P, 'Q', 0, 'water'), rel=1e-9)

    def test_compressor_fed_from_superheated_suction(self):
        """Downstream symptom: compressor outlet must be superheated."""
        from ThermoSim.turbomachinery import Compressor

        M = ThermodynamicModel()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            M.add_point(PENTANE, 'in', T=T10 + 0.022, P=P10, Mass_flowrate=1)
        P_out = 3 * P10
        M.add_point(PENTANE, 'out', P=P_out)
        Compressor(M, 'Comp', 'in', 'out', n_isen=0.8, Calculate=True)
        out = M.Point['out']
        T_sat_out = CP.PropsSI('T', 'P', P_out, 'Q', 1, PENTANE)
        assert out.T > T_sat_out


class TestMixtureNearBubblePoint:
    """For a zeotropic mixture, a T just above the bubble point is a
    two-phase state with tiny quality, NOT vapour: it must stay Q=0."""

    MIX = 'HEOS::R32[0.5]&R125[0.5]'      # R410A, ~0.3 K glide at 10 bar
    P = 1e6

    def setup_method(self):
        config.set_dead_state(T0=298.15, P0=101325)

    def test_just_above_bubble_stays_liquid_side(self):
        T_bub = CP.PropsSI('T', 'P', self.P, 'Q', 0, self.MIX)
        T_dew = CP.PropsSI('T', 'P', self.P, 'Q', 1, self.MIX)
        assert T_dew - T_bub > 0.1          # real glide, wider than the band
        pt, msgs = _make(self.MIX, T_bub + 0.02, self.P)
        assert msgs and 'saturated liquid (Q=0)' in msgs[0]
        assert pt.Q == 0.0
