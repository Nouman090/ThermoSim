"""
test_plotting.py
-----------------
CyclePlotter had no test coverage before this file.

Run with:   pytest tests/test_plotting.py -v
"""

import matplotlib
matplotlib.use('Agg')

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel
from ThermoSim.plotting import CyclePlotter


def _single_loop_model():
    M = ThermodynamicModel()
    M.set_dead_state()
    M.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
    M.add_point('water', '2', P=0.008e6, S=M.Point['1'].S)
    M.add_point('water', '3', P=0.008e6, Q=0)
    M.add_point('water', '4', P=8e6, H=M.Point['3'].H + (8e6 - 0.008e6) / 1000)
    return M


def _cascade_model():
    """Two-loop cascade: R23 low-temp stage, R404a high-temp stage,
    thermally coupled but with distinct fluids and distinct points."""
    M = ThermodynamicModel()
    M.set_dead_state()
    M.add_point('R23', '8', P=1.5e6, T=250, Mass_flowrate=0.5)
    M.add_point('R23', '9', P=1.5e6, Q=1)
    M.add_point('R23', '1', P=0.3e6, T=310)
    M.add_point('R23', '2', P=0.3e6, Q=1)
    M.add_point('R23', '3a', P=0.3e6, T=225)
    M.add_point('R23', '4', P=1.5e6, T=310)
    M.add_point('R404a', '3b', P=0.5e6, T=260)
    M.add_point('R404a', '5', P=1.8e6, T=310)
    M.add_point('R404a', '6', P=1.8e6, Q=1)
    M.add_point('R404a', '7', P=0.5e6, T=260)
    return M


def _dual_pressure_orc_model():
    """Two independent loops, SAME working fluid, no shared points --
    e.g. a dual-pressure ORC's high- and low-pressure stages."""
    M = ThermodynamicModel()
    M.set_dead_state()
    M.add_point('R245fa', 'h1', P=1.4e6, T=380, Mass_flowrate=1.0)
    M.add_point('R245fa', 'h2', P=0.2e6, Q=1)
    M.add_point('R245fa', 'h3', P=0.2e6, Q=0)
    M.add_point('R245fa', 'h4', P=1.4e6, T=310)
    M.add_point('R245fa', 'l1', P=0.8e6, T=350, Mass_flowrate=0.6)
    M.add_point('R245fa', 'l2', P=0.2e6, Q=1)
    M.add_point('R245fa', 'l3', P=0.2e6, Q=0)
    M.add_point('R245fa', 'l4', P=0.8e6, T=305)
    return M


class TestSingleLoopBackwardCompat:
    """A flat point-name list must behave exactly as it did before
    multi-loop support was added: one 'Cycle'-labeled line, one dome."""

    def test_flat_list_gives_single_cycle_label(self):
        M = _single_loop_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(['1', '2', '3', '4'])
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert labels == ['Saturation', 'Cycle']

    def test_flat_list_ts_and_hs_also_unaffected(self):
        M = _single_loop_model()
        plotter = CyclePlotter(M)
        _, ax_ts = plotter.plot_Ts_diagram(['1', '2', '3', '4'])
        assert [t.get_text() for t in ax_ts.get_legend().get_texts()] == \
            ['Saturation', 'Cycle']
        _, ax_hs = plotter.plot_hs_diagram(['1', '2', '3', '4'])
        assert [t.get_text() for t in ax_hs.get_legend().get_texts()] == \
            ['Saturation', 'Cycle']


class TestCascadeMultiLoop:
    """A list of loops -- e.g. the two stages of a cascade system --
    should plot as separate, distinctly colored, distinctly labeled
    cycles with one saturation dome per distinct fluid."""

    def test_two_fluids_two_domes_two_labeled_loops(self):
        M = _cascade_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            [['8', '9', '1', '2', '3a', '4'], ['3b', '5', '6', '7']])
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert 'Saturation (R23)' in labels
        assert 'Saturation (R404a)' in labels
        assert 'Loop 1 (R23)' in labels
        assert 'Loop 2 (R404a)' in labels

    def test_loops_get_distinct_colors(self):
        M = _cascade_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            [['8', '9', '1', '2', '3a', '4'], ['3b', '5', '6', '7']])
        loop_lines = [l for l in ax.get_lines()
                     if (l.get_label() or '').startswith('Loop')]
        assert len(loop_lines) == 2
        assert loop_lines[0].get_color() != loop_lines[1].get_color()

    def test_csv_data_tags_each_row_with_its_loop(self, tmp_path):
        M = _cascade_model()
        plotter = CyclePlotter(M)
        csv_path = str(tmp_path / 'cascade.csv')
        plotter.plot_Ph_diagram(
            [['8', '9', '1', '2', '3a', '4'], ['3b', '5', '6', '7']],
            save_csv=True, csv_path=csv_path)
        import pandas as pd
        df = pd.read_csv(csv_path)
        assert set(df['Loop']) == {1, 2}
        assert set(df[df['Loop'] == 1]['Point']) == {'8', '9', '1', '2', '3a', '4'}
        assert set(df[df['Loop'] == 2]['Point']) == {'3b', '5', '6', '7'}

    def test_shared_point_same_fluid_gives_one_dome(self):
        """Two stages of the SAME fluid sharing a literal point name (e.g.
        a flash-tank / economizer cycle) should still get two loop labels
        but only one deduplicated dome."""
        M = ThermodynamicModel()
        M.set_dead_state()
        M.add_point('R404a', '8', P=1.5e6, T=280)
        M.add_point('R404a', '9', P=1.5e6, Q=1)
        M.add_point('R404a', '1', P=0.3e6, T=310)
        M.add_point('R404a', '2', P=0.3e6, Q=1)
        M.add_point('R404a', '3', P=0.3e6, T=225)
        M.add_point('R404a', '4', P=1.5e6, T=310)
        M.add_point('R404a', '5', P=1.8e6, T=310)
        M.add_point('R404a', '6', P=1.8e6, Q=1)
        M.add_point('R404a', '7', P=0.5e6, T=260)
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            [['8', '9', '1', '2', '3', '4'], ['3', '5', '6', '7']])
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert labels.count('Saturation') == 1
        assert 'Loop 1 (R404a)' in labels
        assert 'Loop 2 (R404a)' in labels


class TestSameFluidLoopLabels:
    """Same-fluid loops (e.g. a dual-pressure ORC) get auto labels that
    only differ by number unless the caller names them explicitly."""

    def test_auto_labels_differ_only_by_number(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram([['h1', 'h2', 'h3', 'h4'],
                                            ['l1', 'l2', 'l3', 'l4']])
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert 'Loop 1 (R245fa)' in labels
        assert 'Loop 2 (R245fa)' in labels

    def test_explicit_loop_labels_override_auto(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            [['h1', 'h2', 'h3', 'h4'], ['l1', 'l2', 'l3', 'l4']],
            loop_labels=['High-Pressure Loop', 'Low-Pressure Loop'])
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert 'High-Pressure Loop' in labels
        assert 'Low-Pressure Loop' in labels
        assert not any('R245fa' in l for l in labels)

    def test_partial_loop_labels_fall_back_to_auto(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            [['h1', 'h2', 'h3', 'h4'], ['l1', 'l2', 'l3', 'l4']],
            loop_labels=['High-Pressure Loop', None])
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert 'High-Pressure Loop' in labels
        assert 'Loop 2 (R245fa)' in labels

    def test_loop_labels_work_on_single_loop_too(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            ['h1', 'h2', 'h3', 'h4'], loop_labels=['My Cycle'])
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert 'My Cycle' in labels
        assert 'Cycle' not in labels


class TestCustomLoopColors:
    """loop_colors lets the caller override the auto palette per loop."""

    def test_explicit_colors_apply(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            [['h1', 'h2', 'h3', 'h4'], ['l1', 'l2', 'l3', 'l4']],
            loop_labels=['High-Pressure Loop', 'Low-Pressure Loop'],
            loop_colors=['crimson', 'navy'])
        handles, labels = ax.get_legend_handles_labels()
        colors = dict(zip(labels, [h.get_color() for h in handles]))
        assert colors['High-Pressure Loop'] == 'crimson'
        assert colors['Low-Pressure Loop'] == 'navy'

    def test_partial_colors_fall_back_to_palette(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            [['h1', 'h2', 'h3', 'h4'], ['l1', 'l2', 'l3', 'l4']],
            loop_colors=[None, 'gold'])
        handles, labels = ax.get_legend_handles_labels()
        colors = dict(zip(labels, [h.get_color() for h in handles]))
        assert colors['Loop 2 (R245fa)'] == 'gold'
        # Loop 1 falls back to the P-h diagram's default first palette color.
        assert colors['Loop 1 (R245fa)'] == 'red'

    def test_single_loop_default_color_unchanged_without_override(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(['h1', 'h2', 'h3', 'h4'])
        handles, labels = ax.get_legend_handles_labels()
        colors = dict(zip(labels, [h.get_color() for h in handles]))
        assert colors['Cycle'] == 'red'

    def test_loop_colors_work_on_single_loop_too(self):
        M = _dual_pressure_orc_model()
        plotter = CyclePlotter(M)
        fig, ax = plotter.plot_Ph_diagram(
            ['h1', 'h2', 'h3', 'h4'], loop_colors=['teal'])
        handles, labels = ax.get_legend_handles_labels()
        colors = dict(zip(labels, [h.get_color() for h in handles]))
        assert colors['Cycle'] == 'teal'
