"""
test_nh3h2o_heat_exchanger.py
-------------------------------
Tests for mixture support in the generic pinch-based HeatExchanger
(HEX_type in {'double_pipe', 'Condenser', 'Evaporator'} -- as opposed to
'SimpleHEX', which never touches this machinery at all).

These exercise the three helpers that make real pinch detection possible
for an NH3-H2O stream (heat_exchangers.py):
  * ``_safe_T_min`` / ``_safe_T_max`` -- root-finding brackets for the
    2-unknown solver, used when the MIXTURE side itself is unknown
  * ``_phase_boundary_alphas``        -- glide-aware breakpoints in the
    pinch T-profile, used by ``_build_profile`` on every solve

Run with:   pytest tests/test_nh3h2o_heat_exchanger.py -v
"""
import sys
import os
import warnings

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, NH3H2OFluid, HeatExchanger

Y_R = 0.995
P_HI = 15e5
P_LO = 4e5


@pytest.fixture(autouse=True)
def _quiet():
    warnings.filterwarnings('ignore')


class TestZeroAndOneUnknown:
    """Mixture side fully known -- exercises _phase_boundary_alphas (glide
    detection) inside _build_profile, called on every solve regardless of
    unknown count."""

    def test_condenser_n1_mixture_hot_side(self):
        M = ThermodynamicModel()
        M.add_point(NH3H2OFluid(Y_R), 'rv', P=P_HI, Q=1.0, Mass_flowrate=0.274)
        M.add_point(NH3H2OFluid(Y_R), 'co', P=P_HI, Q=0.0, Mass_flowrate=0.274)
        M.add_point('water', 'wi', P=2e5, T=298.15)
        M.add_point('water', 'wo', P=2e5, T=305.0)
        hx = HeatExchanger(M, 'CondA', PPT=5.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='rv', Hot_Out_state='co',
                           Cold_In_state='wi', Cold_Out_state='wo',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.Cold_Mass_flowrate > 0
        assert hx.min_dT == pytest.approx(hx.PPT, abs=1e-2) or hx.min_dT > hx.PPT
        assert hx.min_dT >= hx.PPT - 1e-6


class TestTwoUnknown:
    """Mixture side itself unknown -- exercises _safe_T_min/_safe_T_max's
    NH3H2OFluid branch inside the 2-unknown solver's bracket construction."""

    def test_condenser_n2_mixture_Hho_mh(self):
        M = ThermodynamicModel()
        M.add_point(NH3H2OFluid(Y_R), 'hi', P=P_HI, Q=1.0)   # H known
        M.add_point(NH3H2OFluid(Y_R), 'ho', P=P_HI)          # H unknown
        M.add_point('water', 'wi', P=2e5, T=288.15, Mass_flowrate=1.0)
        M.add_point('water', 'wo', P=2e5, T=300.0, Mass_flowrate=1.0)
        hx = HeatExchanger(M, 'CondB', PPT=5.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='wi', Cold_Out_state='wo',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.Hot_Mass_flowrate > 0
        assert hx.min_dT >= hx.PPT - 1e-3
        assert hx.Hot_Out.T < hx.Hot_In.T   # hot stream must cool down


class TestThreeUnknown:
    """Anchor-based 3-unknown solve with the mixture leg fully known and
    the utility leg triply unknown -- exercises the anchor path plus
    repeated _phase_boundary_alphas calls under iterative solving."""

    def test_evaporator_n3_mixture_cold_side(self):
        M = ThermodynamicModel()
        M.add_point(NH3H2OFluid(Y_R), 'ci', P=P_LO, Q=0.10, Mass_flowrate=0.274)
        M.add_point(NH3H2OFluid(Y_R), 'co', P=P_LO)          # H unknown (anchored)
        M.add_point('water', 'hi', P=2e5, T=300.0)           # known
        M.add_point('water', 'ho', P=2e5)                    # H, m unknown
        hx = HeatExchanger(M, 'EvapD', PPT=5.0, HEX_type='Evaporator',
                           HeatAdded=True,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='ci', Cold_Out_state='co',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.Hot_Mass_flowrate > 0
        assert hx.min_dT >= hx.PPT - 1e-3
        assert hx.Hot_Out.T < hx.Hot_In.T


class TestNonRegression:
    """The same 2-unknown configuration with a pure CoolProp fluid must
    behave identically before and after the _safe_T_max refactor (the
    CoolProp branch is unchanged -- CP.PropsSI('TMAX', fluid) verbatim)."""

    def test_pure_fluid_condenser_still_solves(self):
        M = ThermodynamicModel()
        M.add_point('Ammonia', 'hi', P=P_HI, Q=1.0)
        M.add_point('Ammonia', 'ho', P=P_HI)
        M.add_point('water', 'wi', P=2e5, T=288.15, Mass_flowrate=1.0)
        M.add_point('water', 'wo', P=2e5, T=300.0, Mass_flowrate=1.0)
        hx = HeatExchanger(M, 'CondC', PPT=5.0, HEX_type='Condenser',
                           HeatAdded=False,
                           Hot_In_state='hi', Hot_Out_state='ho',
                           Cold_In_state='wi', Cold_Out_state='wo',
                           Calculate=True)
        assert hx.Solution_Status
        assert hx.Hot_Mass_flowrate > 0
