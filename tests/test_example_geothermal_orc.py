"""The executable example of Section 2.3 of the article runs and gives its documented output."""
import os, sys, runpy, warnings

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def test_geothermal_orc_example_runs_and_matches_documented_output():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ns = runpy.run_path(os.path.join(ROOT, 'examples', 'geothermal_orc.py'))
    assert abs(ns['m_wf'] - 4.735) < 0.005          # kg/s
    assert abs(ns['T_b2'] - 104.78) < 0.02          # C
    assert abs(ns['dT_min'] - 8.0) < 0.01           # K
    assert abs(ns['w_net'] - 119.4) < 0.1           # kW
