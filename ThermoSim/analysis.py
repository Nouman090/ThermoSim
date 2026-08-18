"""
analysis.py
-----------
Sensitivity and parametric analysis tools.

Usage:
    from ThermoSim.analysis import SensitivityAnalyzer
    sa = SensitivityAnalyzer(build_function, base_params)
    sa.single_sweep(...)
    sa.multi_sweep(...)
    
    # Or analyze existing CSV data
    sa.analyze_csv('results.csv', ...)
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


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

    def __init__(self, build_func=None, base_params=None):
        self.build_func = build_func
        self.base_params = dict(base_params) if base_params else {}

    # ============================================================== #
    #  CSV I/O Methods
    # ============================================================== #
    def save_results(self, df, filename, overwrite=False):
        """
        Save results DataFrame to CSV.

        Parameters
        ----------
        df : pd.DataFrame
            Results to save.
        filename : str or Path
            Output CSV file path.
        overwrite : bool
            If False, append timestamp to avoid overwriting.

        Returns
        -------
        Path
            Path to saved file.
        """
        filepath = Path(filename)
        
        if not overwrite and filepath.exists():
            # Add timestamp to filename
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            stem = filepath.stem
            suffix = filepath.suffix
            filepath = filepath.parent / f"{stem}_{timestamp}{suffix}"
            warnings.warn(f"File exists. Saving as: {filepath}")
        
        df.to_csv(filepath, index=False)
        print(f"Results saved to: {filepath}")
        return filepath

    def load_results(self, filename):
        """
        Load results from CSV file.

        Parameters
        ----------
        filename : str or Path
            Input CSV file path.

        Returns
        -------
        pd.DataFrame
            Loaded data.
        """
        filepath = Path(filename)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        df = pd.read_csv(filepath)
        print(f"Loaded {len(df)} rows from: {filepath}")
        return df

    def analyze_csv(self, filename, param_name, output_names,
                    x_label=None, figsize=(10, 6), plot=True):
        """
        Analyze existing CSV data (similar to single_sweep but from file).

        Parameters
        ----------
        filename : str or Path
            CSV file containing results.
        param_name : str
            Column name to use as x-axis (parameter varied).
        output_names : list of str
            Column names to plot as outputs.
        x_label : str, optional
            Label for x-axis.
        figsize : tuple
            Figure size.
        plot : bool
            Whether to show plots.

        Returns
        -------
        pd.DataFrame
            Loaded and filtered data.
        """
        df = self.load_results(filename)
        
        # Validate columns
        missing_cols = [col for col in [param_name] + output_names 
                       if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Columns not found in CSV: {missing_cols}")
        
        # Sort by parameter for cleaner plots
        df = df.sort_values(param_name).reset_index(drop=True)
        
        if plot:
            self._plot_single_sweep(df, param_name, output_names, 
                                   x_label or param_name, figsize)
        
        return df

    def analyze_csv_double(self, filename, param1_name, param2_name,
                          output_name, figsize=(10, 8)):
        """
        Analyze 2D sweep data from CSV and create contour plot.

        Parameters
        ----------
        filename : str or Path
            CSV file with 2-parameter sweep results.
        param1_name, param2_name : str
            Column names for the two parameters.
        output_name : str
            Column name for the output to visualize.
        figsize : tuple
            Figure size.

        Returns
        -------
        pd.DataFrame
            Pivot table of results.
        """
        df = self.load_results(filename)
        
        # Validate columns
        required_cols = [param1_name, param2_name, output_name]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Columns not found in CSV: {missing_cols}")
        
        # Create pivot table
        pivot = df.pivot_table(index=param2_name,
                              columns=param1_name,
                              values=output_name)
        
        # Plot
        self._plot_contour(pivot, param1_name, param2_name, 
                          output_name, figsize)
        
        return pivot

    # ============================================================== #
    #  1.  Single-parameter sweep
    # ============================================================== #
    def single_sweep(self, param_name, values, outputs,
                     x_label=None, figsize=(10, 6), plot=True,
                     save_to=None):
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
        figsize : tuple
            Figure size.
        plot : bool
            Whether to show a plot.
        save_to : str or Path, optional
            If provided, save results to this CSV file.

        Returns
        -------
        pd.DataFrame
            Table of results.
        """
        if self.build_func is None:
            raise ValueError("build_func is required for sweep analysis")
        
        x_label = x_label or param_name
        results = {name: [] for name in outputs}
        results[param_name] = []

        for val in values:
            params = dict(self.base_params)
            params[param_name] = val

            try:
                model = self.build_func(params)
                model.ModelSummary(verbose=False)
                results[param_name].append(val)
                for name, func in outputs.items():
                    try:
                        result = func(model)
                        results[name].append(result)
                    except Exception as e:
                        warnings.warn(f"Output '{name}' extraction failed at "
                                    f"{param_name}={val}: {e}")
                        results[name].append(np.nan)
                
                

                
            except Exception as e:
                warnings.warn(f"{param_name}={val} failed: {e}")
                results[param_name].append(val)
                for name in outputs:
                    results[name].append(np.nan)

        df = pd.DataFrame(results)
        
        # Save if requested
        if save_to:
            self.save_results(df, save_to)

        if plot:
            output_names = list(outputs.keys())
            self._plot_single_sweep(df, param_name, output_names, 
                                   x_label, figsize)

        return df

    def _plot_single_sweep(self, df, param_name, output_names, 
                          x_label, figsize):
        """Internal method to plot single sweep results."""
        n_outputs = len(output_names)
        fig, axes = plt.subplots(1, n_outputs,
                                figsize=(figsize[0], figsize[1]),
                                squeeze=False)

        for i, name in enumerate(output_names):
            ax = axes[0][i]
            # Filter out NaN values for plotting
            valid_mask = ~df[name].isna()
            ax.plot(df.loc[valid_mask, param_name], 
                   df.loc[valid_mask, name], 
                   'bo-', linewidth=2, markersize=6)
            ax.set_xlabel(x_label, fontsize=12)
            ax.set_ylabel(name, fontsize=12)
            ax.set_title(f'{name} vs {x_label}', fontsize=13)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    # ============================================================== #
    #  2.  Two-parameter sweep  (contour / heatmap)
    # ============================================================== #
    def double_sweep(self, param1_name, param1_values,
                     param2_name, param2_values,
                     output_name, output_func,
                     figsize=(10, 8), save_to=None):
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
        figsize : tuple
            Figure size.
        save_to : str or Path, optional
            If provided, save results to this CSV file.

        Returns
        -------
        pd.DataFrame
            Full results table (not pivoted).
        """
        if self.build_func is None:
            raise ValueError("build_func is required for sweep analysis")
        
        records = []

        for v1 in param1_values:
            for v2 in param2_values:
                params = dict(self.base_params)
                params[param1_name] = v1
                params[param2_name] = v2

                try:
                    model = self.build_func(params)
                    model.ModelSummary(verbose=False)
                    z = output_func(model)
                except Exception as e:
                    warnings.warn(f"Failed at {param1_name}={v1}, "
                                f"{param2_name}={v2}: {e}")
                    z = np.nan

                records.append({
                    param1_name: v1,
                    param2_name: v2,
                    output_name: z,
                })

        df = pd.DataFrame(records)
        
        # Save if requested
        if save_to:
            self.save_results(df, save_to)

        # Create pivot for plotting
        pivot = df.pivot_table(index=param2_name,
                              columns=param1_name,
                              values=output_name)

        self._plot_contour(pivot, param1_name, param2_name, 
                          output_name, figsize)

        return df

    def _plot_contour(self, pivot, param1_name, param2_name, 
                     output_name, figsize):
        """Internal method to plot contour."""
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get values from pivot table
        param1_values = pivot.columns.values
        param2_values = pivot.index.values
        X, Y = np.meshgrid(param1_values, param2_values)
        Z = pivot.values

        # Check for all NaN
        if np.all(np.isnan(Z)):
            warnings.warn("All values are NaN, cannot create contour plot")
            plt.close(fig)
            return

        # Create contour plot
        contour = ax.contourf(X, Y, Z, levels=20, cmap='viridis')
        fig.colorbar(contour, ax=ax, label=output_name)
        ax.set_xlabel(param1_name, fontsize=13)
        ax.set_ylabel(param2_name, fontsize=13)
        ax.set_title(f'{output_name}', fontsize=15)
        plt.tight_layout()
        plt.show()

    # ============================================================== #
    #  3.  Multi-output sweep (single parameter, many outputs)
    # ============================================================== #
    def multi_output_sweep(self, param_name, values, outputs,
                           x_label=None, figsize=(12, 8), save_to=None):
        """
        Like single_sweep but plots all outputs on one chart with multiple y-axes.

        Parameters
        ----------
        param_name : str
            Parameter to vary.
        values : array-like
            Values to try.
        outputs : dict
            Mapping of display_name → callable(model) → float.
        x_label : str, optional
            Label for x-axis.
        figsize : tuple
            Figure size.
        save_to : str or Path, optional
            If provided, save results to this CSV file.

        Returns
        -------
        pd.DataFrame
            Results table.
        """
        # Use single_sweep but don't plot yet
        df = self.single_sweep(param_name, values, outputs,
                              x_label=x_label, plot=False, save_to=save_to)

        # Create multi-axis plot
        fig, ax1 = plt.subplots(figsize=figsize)
        output_names = list(outputs.keys())
        n_outputs = len(output_names)
        
        # Generate distinct colors
        colors = plt.cm.tab10(np.linspace(0, 1, min(n_outputs, 10)))
        if n_outputs > 10:
            warnings.warn("More than 10 outputs; colors will repeat")
            colors = list(colors) * ((n_outputs // 10) + 1)

        axes = [ax1]
        for i, name in enumerate(output_names):
            if i == 0:
                ax = ax1
            else:
                ax = ax1.twinx()
                # Offset extra y-axes
                if i > 1:
                    ax.spines['right'].set_position(('outward', 60 * (i - 1)))
                axes.append(ax)

            # Filter out NaN values
            valid_mask = ~df[name].isna()
            ax.plot(df.loc[valid_mask, param_name], 
                   df.loc[valid_mask, name],
                   color=colors[i], linewidth=2,
                   marker='o', markersize=5, label=name)
            ax.set_ylabel(name, color=colors[i], fontsize=12)
            ax.tick_params(axis='y', labelcolor=colors[i])

        ax1.set_xlabel(x_label or param_name, fontsize=13)
        ax1.set_title('Multi-Output Sensitivity', fontsize=15)

        # Combined legend
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

    # ============================================================== #
    #  4.  Statistical Analysis
    # ============================================================== #
    def compute_sensitivity_indices(self, df, param_name, output_name):
        """
        Compute basic sensitivity metrics from sweep data.

        Parameters
        ----------
        df : pd.DataFrame
            Results from a sweep.
        param_name : str
            Parameter column name.
        output_name : str
            Output column name.

        Returns
        -------
        dict
            Dictionary with sensitivity metrics.
        """
        # Remove NaN values
        valid_df = df[[param_name, output_name]].dropna()
        
        if len(valid_df) < 2:
            warnings.warn("Insufficient valid data for sensitivity analysis")
            return {}
        
        # Normalize to [0, 1].  A zero-width range would divide by zero and
        # poison every downstream metric with inf/NaN.
        p_range = valid_df[param_name].max() - valid_df[param_name].min()
        o_range = valid_df[output_name].max() - valid_df[output_name].min()
        if p_range == 0:
            warnings.warn(
                f"compute_sensitivity_indices: '{param_name}' does not vary; "
                f"sensitivity is undefined."
            )
            return {}
        param_norm = (valid_df[param_name] - valid_df[param_name].min()) / p_range
        output_norm = ((valid_df[output_name] - valid_df[output_name].min()) / o_range
                       if o_range != 0
                       else valid_df[output_name] * 0.0)
        
        # ── Sensitivity index = mean elasticity ─────────────────────────
        # Up to v3.2.2 this was the mean gradient of output_norm vs
        # param_norm, with BOTH axes rescaled to [0, 1].  That rescaling
        # destroys exactly the information the metric is supposed to carry:
        # for any monotonic, near-linear sweep the normalised rise and run
        # are both 1, so the index came out as ~1.0000 for every parameter
        # and could not rank them against each other.
        #
        # The elasticity (relative sensitivity) is the standard
        # dimensionless alternative:
        #
        #     S = (dY / Y) / (dX / X) = (dY/dX) · (X / Y)
        #
        # It is unit-free, so parameters measured in kelvin and pascals are
        # directly comparable, and it does not collapse to a constant.
        x = valid_df[param_name].to_numpy(dtype=float)
        y = valid_df[output_name].to_numpy(dtype=float)

        dx = np.diff(x)
        dy = np.diff(y)
        with np.errstate(divide='ignore', invalid='ignore'):
            gradients = np.where(dx != 0, dy / dx, np.nan)

        mean_gradient = float(np.nanmean(gradients)) if gradients.size else 0.0

        # Midpoints, so the elasticity is evaluated where the gradient is.
        x_mid = (x[:-1] + x[1:]) / 2.0
        y_mid = (y[:-1] + y[1:]) / 2.0
        with np.errstate(divide='ignore', invalid='ignore'):
            elasticity = np.where(y_mid != 0, gradients * x_mid / y_mid, np.nan)

        sensitivity_index = (float(np.nanmean(np.abs(elasticity)))
                             if elasticity.size else 0.0)

        # Correlation coefficient
        correlation = valid_df[param_name].corr(valid_df[output_name])

        # Range of output
        output_range = valid_df[output_name].max() - valid_df[output_name].min()

        return {
            # Dimensionless elasticity: "a 1 % change in the parameter moves
            # the output by S %".  Comparable across parameters.
            'sensitivity_index': sensitivity_index,
            # Raw slope, in output units per parameter unit.
            'mean_gradient': mean_gradient,
            'correlation': correlation,
            'output_range': output_range,
            'output_min': valid_df[output_name].min(),
            'output_max': valid_df[output_name].max(),
            'param_min': valid_df[param_name].min(),
            'param_max': valid_df[param_name].max(),
        }


# ============================================================== #
#  Standalone utility functions
# ============================================================== #
def compare_csv_files(file1, file2, param_name, output_names, 
                     labels=None, figsize=(12, 6)):
    """
    Compare results from two CSV files side by side.

    Parameters
    ----------
    file1, file2 : str or Path
        Paths to CSV files.
    param_name : str
        Parameter column to use as x-axis.
    output_names : list of str
        Output columns to compare.
    labels : list of str, optional
        Labels for the two datasets (defaults to filenames).
    figsize : tuple
        Figure size.
    """
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    
    if labels is None:
        labels = [Path(file1).stem, Path(file2).stem]
    
    n_outputs = len(output_names)
    fig, axes = plt.subplots(1, n_outputs, figsize=figsize, squeeze=False)
    
    for i, output in enumerate(output_names):
        ax = axes[0][i]
        
        # Plot first dataset
        if output in df1.columns and param_name in df1.columns:
            valid1 = ~df1[output].isna()
            ax.plot(df1.loc[valid1, param_name], df1.loc[valid1, output],
                   'bo-', linewidth=2, label=labels[0])
        
        # Plot second dataset
        if output in df2.columns and param_name in df2.columns:
            valid2 = ~df2[output].isna()
            ax.plot(df2.loc[valid2, param_name], df2.loc[valid2, output],
                   'rs-', linewidth=2, label=labels[1])
        
        ax.set_xlabel(param_name, fontsize=12)
        ax.set_ylabel(output, fontsize=12)
        ax.set_title(f'{output} Comparison', fontsize=13)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()