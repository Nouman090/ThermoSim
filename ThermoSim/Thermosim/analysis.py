"""
analysis.py
-----------
Sensitivity and parametric analysis tools.

Usage:
    from thermocycle.analysis import SensitivityAnalyzer
    sa = SensitivityAnalyzer(Model, build_function)
    sa.single_sweep(...)
    sa.multi_sweep(...)
"""

import copy
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class SensitivityAnalyzer:
    """
    Run parametric / sensitivity studies on a thermodynamic model.

    The key idea: you provide a **build function** that takes a dict
    of parameters and returns a fully-solved ThermodynamicModel.
    This avoids mutating a shared model and keeps every run independent.

    Parameters
    ----------
    build_func : callable
        A function with signature:
            build_func(params: dict) -> ThermodynamicModel
        It should create the model, add points, components, solve it,
        and return the model.
    base_params : dict
        Default parameter values, e.g.
        {'T1': 753.15, 'P1': 8e6, 'eta_t': 0.85, ...}
    """

    def __init__(self, build_func, base_params):
        self.build_func  = build_func
        self.base_params = dict(base_params)

    # ============================================================== #
    #  1.  Single-parameter sweep
    # ============================================================== #
    def single_sweep(self, param_name, values, outputs,
                     x_label=None, figsize=(10, 6), plot=True):
        """
        Vary ONE parameter while keeping everything else at base values.

        Parameters
        ----------
        param_name : str
            Key in base_params to vary, e.g. 'T1'.
        values : array-like
            Array of values to try.
        outputs : dict
            Mapping of  display_name → callable(model) → float
            Example:
                {
                    'Efficiency (%)':   lambda m: m.Efficiency,
                    'Net Power (kW)':   lambda m: m.Net_power / 1e3,
                }
        x_label : str, optional
            Label for the x-axis (defaults to param_name).
        plot : bool
            Whether to show a plot.

        Returns
        -------
        pd.DataFrame
            Table of results.
        """
        x_label = x_label or param_name
        results = {name: [] for name in outputs}
        results[param_name] = []

        for val in values:
            params = dict(self.base_params)
            params[param_name] = val

            try:
                model = self.build_func(params)
                model.ModelSummary()
                results[param_name].append(val)
                for name, func in outputs.items():
                    results[name].append(func(model))
            except Exception as e:
                warnings.warn(f"{param_name}={val} failed: {e}")
                results[param_name].append(val)
                for name in outputs:
                    results[name].append(np.nan)

        df = pd.DataFrame(results)

        if plot:
            n_outputs = len(outputs)
            fig, axes = plt.subplots(1, n_outputs,
                                     figsize=(figsize[0], figsize[1]),
                                     squeeze=False)

            for i, (name, _) in enumerate(outputs.items()):
                ax = axes[0][i]
                ax.plot(df[param_name], df[name], 'bo-',
                        linewidth=2, markersize=6)
                ax.set_xlabel(x_label, fontsize=12)
                ax.set_ylabel(name, fontsize=12)
                ax.set_title(f'{name} vs {x_label}', fontsize=13)
                ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.show()

        return df

    # ============================================================== #
    #  2.  Two-parameter sweep  (contour / heatmap)
    # ============================================================== #
    def double_sweep(self, param1_name, param1_values,
                     param2_name, param2_values,
                     output_name, output_func,
                     figsize=(10, 8)):
        """
        Vary TWO parameters and plot a contour / heatmap.

        Parameters
        ----------
        param1_name, param2_name : str
            Keys in base_params.
        param1_values, param2_values : array-like
            Values to try.
        output_name : str
            Display name for the output.
        output_func : callable(model) -> float
            Function that extracts the output from a solved model.

        Returns
        -------
        pd.DataFrame
            Pivot table of results.
        """
        records = []

        for v1 in param1_values:
            for v2 in param2_values:
                params = dict(self.base_params)
                params[param1_name] = v1
                params[param2_name] = v2

                try:
                    model = self.build_func(params)
                    model.ModelSummary()
                    z = output_func(model)
                except Exception:
                    z = np.nan

                records.append({
                    param1_name: v1,
                    param2_name: v2,
                    output_name: z,
                })

        df = pd.DataFrame(records)
        pivot = df.pivot_table(index=param2_name,
                                columns=param1_name,
                                values=output_name)

        fig, ax = plt.subplots(figsize=figsize)
        X, Y = np.meshgrid(param1_values, param2_values)
        Z = pivot.values

        contour = ax.contourf(X, Y, Z, levels=20, cmap='viridis')
        fig.colorbar(contour, ax=ax, label=output_name)
        ax.set_xlabel(param1_name, fontsize=13)
        ax.set_ylabel(param2_name, fontsize=13)
        ax.set_title(f'{output_name}', fontsize=15)
        plt.tight_layout()
        plt.show()

        return df

    # ============================================================== #
    #  3.  Multi-output sweep (single parameter, many outputs)
    # ============================================================== #
    def multi_output_sweep(self, param_name, values, outputs,
                           x_label=None, figsize=(12, 8)):
        """
        Like single_sweep but plots all outputs on one normalised chart.

        Parameters
        ----------
        Same as single_sweep.

        Returns
        -------
        pd.DataFrame
        """
        df = self.single_sweep(param_name, values, outputs,
                               x_label=x_label, plot=False)

        fig, ax1 = plt.subplots(figsize=figsize)
        colors = plt.cm.tab10(np.linspace(0, 1, len(outputs)))

        axes = [ax1]
        for i, name in enumerate(outputs):
            if i == 0:
                ax = ax1
            else:
                ax = ax1.twinx()
                # offset extra y-axes
                if i > 1:
                    ax.spines['right'].set_position(('outward', 60 * (i - 1)))
                axes.append(ax)

            ax.plot(df[param_name], df[name],
                    color=colors[i], linewidth=2,
                    marker='o', markersize=5, label=name)
            ax.set_ylabel(name, color=colors[i], fontsize=12)
            ax.tick_params(axis='y', labelcolor=colors[i])

        ax1.set_xlabel(x_label or param_name, fontsize=13)
        ax1.set_title('Multi-Output Sensitivity', fontsize=15)

        # combined legend
        lines = []
        labels = []
        for ax in axes:
            ln, lb = ax.get_legend_handles_labels()
            lines.extend(ln)
            labels.extend(lb)
        ax1.legend(lines, labels, loc='best', fontsize=10)

        ax1.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

        return df