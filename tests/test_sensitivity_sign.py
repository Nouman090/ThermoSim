"""
Sign convention of the sensitivity measures (v3.5.0).

``sensitivity_index`` is a magnitude (mean absolute elasticity);
``mean_elasticity`` keeps the sign.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ThermoSim.analysis import SensitivityAnalyzer


def _idx(x, y):
    df = pd.DataFrame({'x': x, 'y': y})
    return SensitivityAnalyzer().compute_sensitivity_indices(df, 'x', 'y')


def test_power_law_elasticity_equals_exponent():
    x = np.linspace(1.0, 2.0, 201)
    r = _idx(x, x ** 2)
    assert abs(r['mean_elasticity'] - 2.0) < 1e-3
    assert abs(r['sensitivity_index'] - 2.0) < 1e-3


def test_decreasing_output_gives_negative_signed_elasticity():
    x = np.linspace(4.0, 12.0, 81)
    y = 200.0 - 10.0 * x                    # falls as x rises
    r = _idx(x, y)
    assert r['mean_elasticity'] < 0
    assert r['sensitivity_index'] > 0       # the index is a magnitude
    assert abs(r['sensitivity_index'] + r['mean_elasticity']) < 1e-12
