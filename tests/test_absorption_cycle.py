"""
test_absorption_cycle.py
-------------------------
Runs validation/validate_nh3h2o_chiller.py -- the single-effect NH3-H2O
absorption chiller built from ThermodynamicModel/Prop/Pump/Expansion_valve/
Mixer/Separator/HeatExchanger -- and relies on its own built-in tolerance
check (it raises AssertionError if any reported quantity deviates from
NH3H2O/demo_cycle.py's reference by more than 0.05%).

Run with:   pytest tests/test_absorption_cycle.py -v
"""
import runpy
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SCRIPT = os.path.join(_REPO_ROOT, 'validation', 'validate_nh3h2o_chiller.py')


def test_absorption_chiller_matches_demo_cycle():
    """Building the cycle from generic ThermoSim components reproduces
    demo_cycle.py's Q_gen, Q_cond, Q_evap, Q_abs, W_pump and COP."""
    runpy.run_path(_SCRIPT, run_name='__main__')
