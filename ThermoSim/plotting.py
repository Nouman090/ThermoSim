"""
plotting.py
-----------
All visualisation lives here — completely separated from computation.

Usage:
    from ThermoSim.plotting import CyclePlotter
    plotter = CyclePlotter(Model)
    plotter.plot_Ts_diagram(['1','2','3','4','1'], save_csv=True, csv_path='Ts_data.csv')
    plotter.plot_Ph_diagram(['1','2','3','4','1'], save_csv=True)
    plotter.plot_hex_profile(Model.Component['HEX1'], save_csv=True)
    plotter.plot_exergy_bar(save_csv=True)
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import warnings
import pandas as pd
import os
from datetime import datetime
from dataclasses import dataclass


@dataclass
class PlotConfig:
    """Centralized plot styling and configuration"""
    cycle_color: str = 'blue'
    cycle_linewidth: float = 2.0
    marker_size: int = 8
    label_offset: tuple = (8, 8)
    grid_alpha: float = 0.3
    
    # Incompressible fluids that don't have saturation domes.
    # NOTE: plain 'Water' was listed here up to v3.2.1.  It is a perfectly
    # ordinary compressible CoolProp fluid (and CoolProp fluid names are
    # case-insensitive), so anyone writing 'Water' instead of 'water' got a
    # silently dome-less diagram.  Only the INCOMP:: backend entry belongs.
    INCOMPRESSIBLE_FLUIDS = {'Therminol66', 'T66', 'DowQ',
                             'INCOMP::T66', 'INCOMP::Water'}
    
    @staticmethod
    def entropy_to_plot(S):
        """Convert entropy to kJ/(kg·K)"""
        return S / 1e3 if S is not None else None
    
    @staticmethod
    def temp_to_celsius(T):
        """Convert temperature to Celsius"""
        return T - 273.15 if T is not None else None
    
    @staticmethod
    def enthalpy_to_plot(H):
        """Convert enthalpy to kJ/kg"""
        return H / 1e3 if H is not None else None
    
    @staticmethod
    def pressure_to_bar(P):
        """Convert pressure to bar"""
        return P / 1e5 if P is not None else None
    
    @staticmethod
    def power_to_kw(W):
        """Convert power to kW"""
        return W / 1e3 if W is not None else None


class CyclePlotter:
    """Handles every kind of plot for a ThermodynamicModel."""

    def __init__(self, model):
        """
        Parameters
        ----------
        model : ThermodynamicModel
            A solved model whose Point and Component dicts are populated.
        """
        if not hasattr(model, 'Point') or not hasattr(model, 'Component'):
            raise ValueError("Model must have 'Point' and 'Component' attributes")
        
        self.model = model
        self.config = PlotConfig()

    # ================================================================ #
    #  Utility Methods
    # ================================================================ #
    def _safe_get_property(self, point, prop_name, default=None):
        """Safely retrieve property with fallback"""
        try:
            value = getattr(point, prop_name)
            return value if value is not None else default
        except AttributeError:
            warnings.warn(f"Point has no attribute '{prop_name}'")
            return default

    def _can_draw_dome(self, fluid):
        """Check if saturation dome can be drawn for this fluid"""
        if not fluid or fluid in self.config.INCOMPRESSIBLE_FLUIDS:
            return False
        try:
            CP.PropsSI('TCRIT', fluid)
            return True
        except Exception:
            return False

    # ================================================================ #
    #  Process paths
    # ================================================================ #
    #  A cycle is not a polygon.  The line joining two state points is the
    #  *process* between them, and most processes are curves:
    #
    #    - Isobaric heating crosses the dome, so it runs
    #        sub-cooled -> saturated liquid -> saturated vapour -> superheated
    #      with a sharp kink at each saturation crossing and a dead-flat
    #      section in between (constant T while boiling).
    #    - Throttling is isenthalpic, not a straight line in T-s.
    #    - Real expansion and compression curve away from the isentrope.
    #
    #  Both endpoints always have a well-defined (P, H), and (P, H) fixes a
    #  state everywhere -- unlike (P, T), which is degenerate inside the
    #  dome.  So we interpolate P and H linearly and let CoolProp supply
    #  T and S at each step.  Constant-P and constant-H processes fall out
    #  of that automatically, and saturation crossings are inserted exactly
    #  so the kinks stay sharp instead of being rounded off by sampling.

    def _saturation_crossings(self, fluid, P_a, H_a, P_b, H_b):
        """
        Fractions along an isobaric segment where H crosses h_f or h_g.

        Returns [] for non-isobaric segments, supercritical pressures, or
        any fluid without a dome.
        """
        if not self._can_draw_dome(fluid):
            return []
        # only meaningful when the pressure is (near enough) constant
        if abs(P_b - P_a) > 1e-6 * max(abs(P_a), abs(P_b), 1.0):
            return []
        if abs(H_b - H_a) < 1e-9:
            return []
        try:
            if P_a >= CP.PropsSI('PCRIT', fluid):
                return []
            h_f = CP.PropsSI('H', 'P', P_a, 'Q', 0, fluid)
            h_g = CP.PropsSI('H', 'P', P_a, 'Q', 1, fluid)
        except Exception:
            return []

        fracs = []
        for h_sat in (h_f, h_g):
            f = (h_sat - H_a) / (H_b - H_a)
            if 1e-9 < f < 1 - 1e-9:
                # duplicate the crossing so the kink is a true corner
                fracs.extend([f - 1e-9, f + 1e-9])
        return fracs

    def _find_component(self, name_a, name_b):
        """The component whose inlet is name_a and outlet is name_b, if any."""
        for comp in self.model.Component.values():
            outs = getattr(comp, 'Out_states', None)
            out_ok = (getattr(comp, 'Out_state', None) == name_b or
                      (outs is not None and name_b in outs))
            ins = getattr(comp, 'In_states', None)
            in_ok = (getattr(comp, 'In_state', None) == name_a or
                     (ins is not None and name_a in ins))
            if in_ok and out_ok:
                return comp
        return None

    def _work_path(self, fluid, P_a, H_a, S_a, P_b, H_b, fracs):
        """
        Expansion / compression along a constant-efficiency condition line.

        Interpolating (P, H) linearly is wrong for turbomachinery: it makes
        entropy *fall* in mid-expansion, which no adiabatic machine does.
        The standard construction instead applies the same isentropic
        efficiency at every intermediate pressure, so entropy rises
        monotonically and the endpoints are reproduced exactly.

        The efficiency is derived from the endpoints rather than read off
        the component, so it stays correct even when the outlet state was
        supplied directly.
        """
        try:
            H_s_out = CP.PropsSI('H', 'P', P_b, 'S', S_a, fluid)
        except Exception:
            return None
        ideal = H_s_out - H_a
        if abs(ideal) < 1e-6:
            return None
        eta = (H_b - H_a) / ideal          # <1 expansion, >1 compression

        P_path, H_path = [], []
        for f in fracs:
            P = P_a + f * (P_b - P_a)
            try:
                H = H_a + eta * (CP.PropsSI('H', 'P', P, 'S', S_a, fluid) - H_a)
            except Exception:
                continue
            P_path.append(P)
            H_path.append(H)
        return (P_path, H_path) if len(P_path) >= 3 else None

    def _process_path(self, pt_a, pt_b, n=60, name_a=None, name_b=None):
        """
        Trace the real process from ``pt_a`` to ``pt_b``.

        Returns a dict of equal-length lists: P (Pa), H (J/kg), T (K),
        S (J/kg/K).  Falls back to the two endpoints alone -- i.e. the old
        straight line -- whenever the path can't be evaluated.
        """
        P_a = self._safe_get_property(pt_a, 'P')
        P_b = self._safe_get_property(pt_b, 'P')
        H_a = self._safe_get_property(pt_a, 'H')
        H_b = self._safe_get_property(pt_b, 'H')
        T_a = self._safe_get_property(pt_a, 'T')
        T_b = self._safe_get_property(pt_b, 'T')
        S_a = self._safe_get_property(pt_a, 'S')
        S_b = self._safe_get_property(pt_b, 'S')

        endpoints = {'P': [P_a, P_b], 'H': [H_a, H_b],
                     'T': [T_a, T_b], 'S': [S_a, S_b]}

        if None in (P_a, P_b, H_a, H_b):
            return endpoints

        fluid = getattr(pt_a, 'fluid', None)
        if fluid is None or fluid != getattr(pt_b, 'fluid', None):
            return endpoints
        # Incompressible correlations are handled by Prop, not CoolProp's
        # PropsSI, and they never cross a dome -- a straight line is right.
        if fluid in self.config.INCOMPRESSIBLE_FLUIDS:
            return endpoints

        fracs = sorted(set(
            list(np.linspace(0.0, 1.0, max(int(n), 2))) +
            self._saturation_crossings(fluid, P_a, H_a, P_b, H_b)
        ))

        # Work-transfer machines follow a condition line, not a chord.
        path = None
        comp = self._find_component(name_a, name_b) if name_a and name_b else None
        is_work = type(comp).__name__ in ('Turbine', 'Pump', 'Compressor')
        if is_work and S_a is not None and abs(P_b - P_a) > 0:
            path = self._work_path(fluid, P_a, H_a, S_a, P_b, H_b, fracs)

        if path is None:
            # Everything else -- isobaric heating, throttling, pipes --
            # is traced by interpolating (P, H), which reproduces constant
            # pressure and constant enthalpy exactly.
            path = ([P_a + f * (P_b - P_a) for f in fracs],
                    [H_a + f * (H_b - H_a) for f in fracs])

        P_path, H_path, T_path, S_path = [], [], [], []
        for P, H in zip(*path):
            try:
                T = CP.PropsSI('T', 'P', P, 'H', H, fluid)
                S = CP.PropsSI('S', 'P', P, 'H', H, fluid)
            except Exception:
                continue
            P_path.append(P)
            H_path.append(H)
            T_path.append(T)
            S_path.append(S)

        # Need a real curve, and the endpoints must be the ones the model
        # solved -- not CoolProp's round-trip of them.
        if len(P_path) < 3:
            return endpoints
        for path_, a, b in ((P_path, P_a, P_b), (H_path, H_a, H_b),
                            (T_path, T_a, T_b), (S_path, S_a, S_b)):
            if a is not None:
                path_[0] = a
            if b is not None:
                path_[-1] = b
        return {'P': P_path, 'H': H_path, 'T': T_path, 'S': S_path}

    def _cycle_path(self, points, x_key, y_key, n=60):
        """
        Concatenate the process paths across a whole loop.

        ``points`` is a list of (name, Prop).  Returns (x, y) in SI.
        """
        xs, ys = [], []
        for (na, pt_a), (nb, pt_b) in zip(points, points[1:]):
            seg = self._process_path(pt_a, pt_b, n=n, name_a=na, name_b=nb)
            x, y = seg[x_key], seg[y_key]
            if None in x or None in y:
                x, y = [x[0], x[-1]], [y[0], y[-1]]
            if xs:                      # avoid duplicating the shared node
                x, y = x[1:], y[1:]
            xs.extend(x)
            ys.extend(y)
        return xs, ys

    def _generate_csv_filename(self, base_name, csv_path=None):
        """Generate CSV filename with timestamp if not provided"""
        if csv_path:
            return csv_path
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"{base_name}_{timestamp}.csv"

    def _save_dataframe_to_csv(self, df, csv_path, plot_type):
        """Save DataFrame to CSV with error handling"""
        try:
            # Create directory if it doesn't exist
            directory = os.path.dirname(csv_path)
            if directory and not os.path.exists(directory):
                os.makedirs(directory)
            
            df.to_csv(csv_path, index=False)
            print(f"✓ {plot_type} data saved to: {csv_path}")
            return True
        except Exception as e:
            warnings.warn(f"Failed to save CSV: {e}")
            return False

    # ================================================================ #
    #  1.  T-s Diagram
    # ================================================================ #
    def plot_Ts_diagram(self, loop_points=None, fluid=None,
                        show_dome=True, show_labels=True,
                        title='T-s Diagram', figsize=(10, 7),
                        save_csv=False, csv_path=None,
                        process_path=True, path_N=60):
        """
        Plot a Temperature–Entropy diagram.

        Parameters
        ----------
        loop_points : list of str, optional
            State-point names in cycle order, e.g. ['1','2','3','4','1'].
            The last name should repeat the first to close the loop.
            If None, every point is plotted as a scatter dot.
        fluid : str, optional
            Fluid name for the saturation dome.  Auto-detected from
            the first point in loop_points when omitted.
        show_dome : bool
            Draw the saturation dome behind the cycle.
        show_labels : bool
            Annotate each state-point number on the diagram.
        title : str
            Plot title.
        figsize : tuple
            Figure size in inches.
        save_csv : bool
            If True, save plot data to CSV file.
        csv_path : str, optional
            Custom CSV file path. If None, auto-generates with timestamp.
        process_path : bool
            Trace the real thermodynamic path between state points instead
            of joining them with straight lines.  An isobaric segment then
            shows its saturation kinks (sub-cooled -> saturated liquid ->
            saturated vapour -> superheated) and a throttle stays
            isenthalpic.  Set False for the old polygon.
        path_N : int
            Points sampled per process segment when ``process_path`` is on.

        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        fig, ax = plt.subplots(figsize=figsize)

        # ---- Auto-detect fluid ----
        if fluid is None and loop_points:
            try:
                fluid = self.model.Point[loop_points[0]].fluid
            except (KeyError, AttributeError):
                warnings.warn("Could not auto-detect fluid")

        # ---- Saturation dome ----
        if show_dome and self._can_draw_dome(fluid):
            self._draw_saturation_dome(ax, fluid)

        # ---- Prepare data collection for CSV ----
        csv_data = []

        # ---- Cycle plotting ----
        if loop_points is not None:
            S_vals, T_vals, valid_names = [], [], []
            
            for name in loop_points:
                try:
                    pt = self.model.Point[name]
                except KeyError:
                    warnings.warn(f"Point '{name}' not found in model")
                    continue
                
                S = self._safe_get_property(pt, 'S')
                T = self._safe_get_property(pt, 'T')
                
                if S is None or T is None:
                    warnings.warn(f"Point '{name}' missing S or T")
                    continue
                
                S_plot = S / 1e3
                T_plot = T - 273.15
                
                S_vals.append(S_plot)
                T_vals.append(T_plot)
                valid_names.append(name)
                
                # Collect data for CSV
                csv_data.append({
                    'Point': name,
                    'Entropy_kJ_per_kg_K': S_plot,
                    'Temperature_C': T_plot,
                    'Entropy_J_per_kg_K': S,
                    'Temperature_K': T
                })
            
            if len(S_vals) < 2:
                warnings.warn("Not enough valid points to draw cycle")
                return fig, ax

            if process_path:
                pts = [(n_, self.model.Point[n_]) for n_ in valid_names]
                S_line, T_line = self._cycle_path(pts, 'S', 'T', n=path_N)
                S_line = [s / 1e3 for s in S_line]
                T_line = [t - 273.15 for t in T_line]
            else:
                S_line, T_line = S_vals, T_vals

            ax.plot(S_line, T_line, '-', color=self.config.cycle_color,
                    linewidth=self.config.cycle_linewidth,
                    zorder=4, label='Cycle')
            ax.plot(S_vals, T_vals, 'o', color=self.config.cycle_color,
                    markersize=self.config.marker_size, zorder=5)
            
            if show_labels:
                # Check if cycle is closed
                is_closed = (valid_names[0] == valid_names[-1])
                label_points = valid_names[:-1] if is_closed else valid_names
                
                for i, name in enumerate(label_points):
                    ax.annotate(
                        name, (S_vals[i], T_vals[i]),
                        textcoords='offset points',
                        xytext=(8, 8), fontsize=12,
                        fontweight='bold', color='darkblue',
                        bbox=dict(boxstyle='round,pad=0.3', 
                                 facecolor='yellow', alpha=0.3)
                    )
        else:
            # Scatter mode
            for name, pt in self.model.Point.items():
                S = self._safe_get_property(pt, 'S')
                T = self._safe_get_property(pt, 'T')
                
                if S is not None and T is not None:
                    S_plot = S / 1e3
                    T_plot = T - 273.15
                    
                    ax.scatter(S_plot, T_plot, s=60, zorder=5)
                    
                    if show_labels:
                        ax.annotate(name, (S_plot, T_plot),
                                  textcoords='offset points',
                                  xytext=(5, 5), fontsize=10)
                    
                    # Collect data for CSV
                    csv_data.append({
                        'Point': name,
                        'Entropy_kJ_per_kg_K': S_plot,
                        'Temperature_C': T_plot,
                        'Entropy_J_per_kg_K': S,
                        'Temperature_K': T
                    })

        ax.set_xlabel('Entropy  [kJ/(kg·K)]', fontsize=13)
        ax.set_ylabel('Temperature  [°C]', fontsize=13)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.3)
        if loop_points:
            ax.legend()
        plt.tight_layout()
        plt.show()

        # ---- Save CSV if requested ----
        if save_csv and csv_data:
            df = pd.DataFrame(csv_data)
            csv_filename = self._generate_csv_filename('Ts_diagram', csv_path)
            self._save_dataframe_to_csv(df, csv_filename, 'T-s diagram')
        
        return fig, ax

    # ================================================================ #
    #  2.  P-h Diagram
    # ================================================================ #
    def plot_Ph_diagram(self, loop_points=None, fluid=None,
                        show_dome=True, show_labels=True,
                        title='P-h Diagram', figsize=(10, 7),
                        save_csv=False, csv_path=None,
                        process_path=True, path_N=60):
        """
        Plot a log-Pressure vs Enthalpy diagram.

        Parameters
        ----------
        loop_points : list of str, optional
            State-point names in cycle order.
        fluid : str, optional
            Fluid name for saturation dome.
        show_dome : bool
            Draw the saturation dome.
        show_labels : bool
            Annotate state points.
        title : str
            Plot title.
        figsize : tuple
            Figure size.
        save_csv : bool
            If True, save plot data to CSV.
        csv_path : str, optional
            Custom CSV file path.

        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        fig, ax = plt.subplots(figsize=figsize)

        if fluid is None and loop_points:
            try:
                fluid = self.model.Point[loop_points[0]].fluid
            except (KeyError, AttributeError):
                warnings.warn("Could not auto-detect fluid")

        # dome
        if show_dome and self._can_draw_dome(fluid):
            self._draw_saturation_dome_Ph(ax, fluid)

        csv_data = []

        if loop_points is not None:
            H_vals, P_vals, valid_names = [], [], []
            
            for name in loop_points:
                try:
                    pt = self.model.Point[name]
                except KeyError:
                    warnings.warn(f"Point '{name}' not found")
                    continue
                
                H = self._safe_get_property(pt, 'H')
                P = self._safe_get_property(pt, 'P')
                
                if H is None or P is None:
                    warnings.warn(f"Point '{name}' missing H or P")
                    continue
                
                H_plot = H / 1e3
                P_plot = P / 1e5
                
                H_vals.append(H_plot)
                P_vals.append(P_plot)
                valid_names.append(name)
                
                csv_data.append({
                    'Point': name,
                    'Enthalpy_kJ_per_kg': H_plot,
                    'Pressure_bar': P_plot,
                    'Enthalpy_J_per_kg': H,
                    'Pressure_Pa': P
                })
            
            if len(H_vals) < 2:
                warnings.warn("Not enough valid points to draw cycle")
                return fig, ax

            if process_path:
                pts = [(n_, self.model.Point[n_]) for n_ in valid_names]
                H_line, P_line = self._cycle_path(pts, 'H', 'P', n=path_N)
                H_line = [h / 1e3 for h in H_line]
                P_line = [p / 1e5 for p in P_line]
            else:
                H_line, P_line = H_vals, P_vals

            ax.semilogy(H_line, P_line, '-', color='red',
                        linewidth=2, zorder=4, label='Cycle')
            ax.semilogy(H_vals, P_vals, 'o', color='red',
                        markersize=8, zorder=5)

            if show_labels:
                is_closed = (valid_names[0] == valid_names[-1])
                label_points = valid_names[:-1] if is_closed else valid_names
                
                for i, name in enumerate(label_points):
                    ax.annotate(
                        name, (H_vals[i], P_vals[i]),
                        textcoords='offset points',
                        xytext=(8, 8), fontsize=12,
                        fontweight='bold', color='darkred',
                        bbox=dict(boxstyle='round,pad=0.3',
                                 facecolor='yellow', alpha=0.3)
                    )
        else:
            for name, pt in self.model.Point.items():
                H = self._safe_get_property(pt, 'H')
                P = self._safe_get_property(pt, 'P')
                
                if H is not None and P is not None:
                    H_plot = H / 1e3
                    P_plot = P / 1e5
                    
                    ax.scatter(H_plot, P_plot, s=60, zorder=5)
                    
                    if show_labels:
                        ax.annotate(name, (H_plot, P_plot),
                                  textcoords='offset points',
                                  xytext=(5, 5), fontsize=10)
                    
                    csv_data.append({
                        'Point': name,
                        'Enthalpy_kJ_per_kg': H_plot,
                        'Pressure_bar': P_plot,
                        'Enthalpy_J_per_kg': H,
                        'Pressure_Pa': P
                    })

        ax.set_xlabel('Enthalpy  [kJ/kg]', fontsize=13)
        ax.set_ylabel('Pressure  [bar]', fontsize=13)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.3, which='both')
        if loop_points:
            ax.legend()
        plt.tight_layout()
        plt.show()

        # Save CSV
        if save_csv and csv_data:
            df = pd.DataFrame(csv_data)
            csv_filename = self._generate_csv_filename('Ph_diagram', csv_path)
            self._save_dataframe_to_csv(df, csv_filename, 'P-h diagram')

        return fig, ax

    # ================================================================ #
    #  3.  h-s Diagram  (Mollier)
    # ================================================================ #
    def plot_hs_diagram(self, loop_points=None, fluid=None,
                        show_dome=True, show_labels=True,
                        title='h-s Diagram', figsize=(10, 7),
                        save_csv=False, csv_path=None,
                        process_path=True, path_N=60):
        """
        Enthalpy vs Entropy (Mollier) diagram.

        Parameters
        ----------
        loop_points : list of str, optional
            State-point names in cycle order.
        fluid : str, optional
            Fluid name for saturation dome.
        show_dome : bool
            Draw saturation dome.
        show_labels : bool
            Annotate points.
        title : str
            Plot title.
        figsize : tuple
            Figure size.
        save_csv : bool
            If True, save plot data to CSV.
        csv_path : str, optional
            Custom CSV file path.

        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        fig, ax = plt.subplots(figsize=figsize)

        if fluid is None and loop_points:
            try:
                fluid = self.model.Point[loop_points[0]].fluid
            except (KeyError, AttributeError):
                warnings.warn("Could not auto-detect fluid")

        if show_dome and self._can_draw_dome(fluid):
            self._draw_saturation_dome_hs(ax, fluid)

        csv_data = []

        if loop_points is not None:
            S_vals, H_vals, valid_names = [], [], []
            
            for n in loop_points:
                try:
                    pt = self.model.Point[n]
                except KeyError:
                    warnings.warn(f"Point '{n}' not found")
                    continue
                
                S = self._safe_get_property(pt, 'S')
                H = self._safe_get_property(pt, 'H')
                
                if S is None or H is None:
                    warnings.warn(f"Point '{n}' missing S or H")
                    continue
                
                S_plot = S / 1e3
                H_plot = H / 1e3
                
                S_vals.append(S_plot)
                H_vals.append(H_plot)
                valid_names.append(n)
                
                csv_data.append({
                    'Point': n,
                    'Entropy_kJ_per_kg_K': S_plot,
                    'Enthalpy_kJ_per_kg': H_plot,
                    'Entropy_J_per_kg_K': S,
                    'Enthalpy_J_per_kg': H
                })
            
            if len(S_vals) < 2:
                warnings.warn("Not enough valid points to draw cycle")
                return fig, ax
            
            if process_path:
                pts = [(n_, self.model.Point[n_]) for n_ in valid_names]
                S_line, H_line = self._cycle_path(pts, 'S', 'H', n=path_N)
                S_line = [s / 1e3 for s in S_line]
                H_line = [h / 1e3 for h in H_line]
            else:
                S_line, H_line = S_vals, H_vals

            ax.plot(S_line, H_line, '-', color='green',
                    linewidth=2, zorder=4, label='Cycle')
            ax.plot(S_vals, H_vals, 'o', color='green',
                    markersize=8, zorder=5)

            if show_labels:
                is_closed = (valid_names[0] == valid_names[-1])
                label_points = valid_names[:-1] if is_closed else valid_names
                
                for i, name in enumerate(label_points):
                    ax.annotate(
                        name, (S_vals[i], H_vals[i]),
                        textcoords='offset points',
                        xytext=(8, 8), fontsize=12,
                        fontweight='bold', color='darkgreen',
                        bbox=dict(boxstyle='round,pad=0.3',
                                 facecolor='yellow', alpha=0.3)
                    )
        else:
            for name, pt in self.model.Point.items():
                S = self._safe_get_property(pt, 'S')
                H = self._safe_get_property(pt, 'H')
                
                if S is not None and H is not None:
                    S_plot = S / 1e3
                    H_plot = H / 1e3
                    
                    ax.scatter(S_plot, H_plot, s=60, zorder=5)
                    
                    if show_labels:
                        ax.annotate(name, (S_plot, H_plot),
                                  textcoords='offset points',
                                  xytext=(5, 5), fontsize=10)
                    
                    csv_data.append({
                        'Point': name,
                        'Entropy_kJ_per_kg_K': S_plot,
                        'Enthalpy_kJ_per_kg': H_plot,
                        'Entropy_J_per_kg_K': S,
                        'Enthalpy_J_per_kg': H
                    })

        ax.set_xlabel('Entropy  [kJ/(kg·K)]', fontsize=13)
        ax.set_ylabel('Enthalpy  [kJ/kg]', fontsize=13)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.3)
        if loop_points:
            ax.legend()
        plt.tight_layout()
        plt.show()

        # Save CSV
        if save_csv and csv_data:
            df = pd.DataFrame(csv_data)
            csv_filename = self._generate_csv_filename('hs_diagram', csv_path)
            self._save_dataframe_to_csv(df, csv_filename, 'h-s diagram')

        return fig, ax

    # ================================================================ #
    #  4.  Heat-Exchanger Temperature Profile
    # ================================================================ #
    def plot_hex_profile(self, hex_id, div_N=200, figsize=(12, 5),
                        save_csv=False, csv_path=None):
        """
        Plot the temperature profile and ΔT along a heat exchanger.

        Parameters
        ----------
        hex_id : str
            The ID of the HeatExchanger component in the model.
        div_N : int
            Number of discretisation segments.
        figsize : tuple
            Figure size.
        save_csv : bool
            If True, save profile data to CSV.
        csv_path : str, optional
            Custom CSV file path.

        Returns
        -------
        min_dT : float
            The minimum temperature difference (pinch) in K.
        """
        try:
            hx = self.model.Component[hex_id]
        except KeyError:
            raise ValueError(f"Component '{hex_id}' not found in model")
        
        if not hx.Solution_Status:
            raise ValueError(f"{hex_id} has not been solved yet.")

        if hx.HEX_type == 'SimpleHEX':
            raise ValueError(
                "SimpleHEX has only two points — no internal profile."
            )

        Q_total = hx.Q

        if Q_total is None:
            raise ValueError(
                f"{hex_id} has no heat duty (Q is None). Solve the component "
                f"before plotting its profile."
            )
        if Q_total == 0:
            warnings.warn(f"{hex_id} has zero heat transfer - profile will be flat")
            return 0.0

        if hx.Hot_Mass_flowrate in (None, 0) or hx.Cold_Mass_flowrate in (None, 0):
            raise ValueError(
                f"{hex_id} has a missing or zero mass flow rate "
                f"(hot={hx.Hot_Mass_flowrate}, cold={hx.Cold_Mass_flowrate})."
            )

        q = Q_total / div_N
        dP_h = hx.Hot_In.P - hx.Hot_Out.P
        dP_c = hx.Cold_In.P - hx.Cold_Out.P

        Th = np.zeros(div_N + 1)
        Tc = np.zeros(div_N + 1)
        Q_axis = np.linspace(0, Q_total / 1e3, div_N + 1)

        csv_data = []

        for n in range(div_N + 1):
            hh = hx.Hot_Out.H + q * n / hx.Hot_Mass_flowrate
            hc = hx.Cold_In.H + q * n / hx.Cold_Mass_flowrate
            Ph = hx.Hot_Out.P + (dP_h / div_N) * n
            Pc = hx.Cold_In.P - (dP_c / div_N) * n

            T_hot = self.model.Prop(
                hx.Hot_In.fluid, StatePointName='_plot_h',
                H=hh, P=Ph
            ).T - 273.15
            
            T_cold = self.model.Prop(
                hx.Cold_In.fluid, StatePointName='_plot_c',
                H=hc, P=Pc
            ).T - 273.15
            
            Th[n] = T_hot
            Tc[n] = T_cold
            
            csv_data.append({
                'Heat_Transferred_kW': Q_axis[n],
                'Hot_Temperature_C': T_hot,
                'Cold_Temperature_C': T_cold,
                'Delta_T_K': T_hot - T_cold,
                'Hot_Enthalpy_J_per_kg': hh,
                'Cold_Enthalpy_J_per_kg': hc,
                'Hot_Pressure_Pa': Ph,
                'Cold_Pressure_Pa': Pc
            })

        dT = Th - Tc
        min_dT = min(dT)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # --- temperature profile ---
        ax1.plot(Q_axis, Th, 'r-', linewidth=2, label='Hot side')
        ax1.plot(Q_axis, Tc, 'b-', linewidth=2, label='Cold side')
        ax1.fill_between(Q_axis, Tc, Th, alpha=0.08, color='orange')
        ax1.set_xlabel('Heat transferred  [kW]', fontsize=12)
        ax1.set_ylabel('Temperature  [°C]', fontsize=12)
        ax1.set_title(f'{hex_id} — Temperature Profile', fontsize=13)
        ax1.legend(fontsize=11)
        ax1.grid(True, alpha=0.3)

        # --- ΔT profile ---
        ax2.plot(Q_axis, dT, 'g-', linewidth=2)
        ax2.axhline(y=min_dT, color='red', linestyle='--',
                     linewidth=1.5,
                     label=f'Min ΔT = {min_dT:.2f} K')
        idx_min = np.argmin(dT)
        ax2.plot(Q_axis[idx_min], dT[idx_min], 'rv', markersize=12)
        ax2.set_xlabel('Heat transferred  [kW]', fontsize=12)
        ax2.set_ylabel('ΔT  [K]', fontsize=12)
        ax2.set_title(f'{hex_id} — Pinch Diagram', fontsize=13)
        ax2.legend(fontsize=11)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        # Save CSV
        if save_csv and csv_data:
            df = pd.DataFrame(csv_data)
            csv_filename = self._generate_csv_filename(f'{hex_id}_profile', csv_path)
            self._save_dataframe_to_csv(df, csv_filename, f'{hex_id} temperature profile')

        return min_dT

    # ================================================================ #
    #  5.  Exergy Destruction Bar Chart
    # ================================================================ #
    def plot_exergy_bar(self, figsize=(12, 6), as_percentage=False,
                       save_csv=False, csv_path=None):
        """
        Bar chart of exergy destruction for every component.

        Parameters
        ----------
        figsize : tuple
            Figure size.
        as_percentage : bool
            If True, show each component's share of total Ex_D.
        save_csv : bool
            If True, save exergy data to CSV.
        csv_path : str, optional
            Custom CSV file path.

        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        names = []
        values = []
        csv_data = []
        
        for comp_id, comp in self.model.Component.items():
            if isinstance(comp.Ex_D, (int, float)):
                names.append(comp_id)
                values.append(comp.Ex_D)

        if not values:
            print("No numeric exergy destruction data found.")
            return None, None

        values = np.array(values)
        total = values.sum()

        fig, ax = plt.subplots(figsize=figsize)
        colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(values)))

        if as_percentage and total > 0:
            pct = values / total * 100
            bars = ax.bar(names, pct, color=colors, edgecolor='black',
                          linewidth=0.5)
            ax.set_ylabel('Exergy Destruction  [%]', fontsize=13)

            for i, (bar, p, name) in enumerate(zip(bars, pct, names)):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.5,
                        f'{p:.1f}%', ha='center', fontsize=10)
                
                csv_data.append({
                    'Component': name,
                    'Exergy_Destruction_W': values[i],
                    'Exergy_Destruction_kW': values[i] / 1e3,
                    'Percentage': p
                })
        else:
            bars = ax.bar(names, values / 1e3, color=colors,
                          edgecolor='black', linewidth=0.5)
            ax.set_ylabel('Exergy Destruction  [kW]', fontsize=13)

            for i, (bar, v, name) in enumerate(zip(bars, values / 1e3, names)):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(values/1e3)*0.01,
                        f'{v:.2f}', ha='center', fontsize=10)
                
                pct = (values[i] / total * 100) if total > 0 else 0
                csv_data.append({
                    'Component': name,
                    'Exergy_Destruction_W': values[i],
                    'Exergy_Destruction_kW': v,
                    'Percentage': pct
                })

        ax.set_title('Component-wise Exergy Destruction', fontsize=15)
        plt.xticks(rotation=45, ha='right', fontsize=11)
        plt.tight_layout()
        plt.show()

        # Save CSV
        if save_csv and csv_data:
            df = pd.DataFrame(csv_data)
            csv_filename = self._generate_csv_filename('exergy_destruction', csv_path)
            self._save_dataframe_to_csv(df, csv_filename, 'Exergy destruction')

        return fig, ax

    # ================================================================ #
    #  6.  Exergy Pie Chart
    # ================================================================ #
    def plot_exergy_pie(self, figsize=(8, 8), save_csv=False, csv_path=None):
        """
        Pie chart showing share of exergy destruction.

        Parameters
        ----------
        figsize : tuple
            Figure size.
        save_csv : bool
            If True, save data to CSV.
        csv_path : str, optional
            Custom CSV file path.

        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        names = []
        values = []
        csv_data = []
        
        for comp_id, comp in self.model.Component.items():
            if isinstance(comp.Ex_D, (int, float)) and comp.Ex_D > 0:
                names.append(comp_id)
                values.append(comp.Ex_D)

        if not values:
            print("No numeric exergy destruction data found.")
            return None, None

        total = sum(values)
        
        fig, ax = plt.subplots(figsize=figsize)
        colors = plt.cm.Set3(np.linspace(0, 1, len(values)))
        wedges, texts, autotexts = ax.pie(
            values, labels=names, autopct='%1.1f%%',
            colors=colors, startangle=90,
            textprops={'fontsize': 11}
        )
        ax.set_title('Exergy Destruction Shares', fontsize=15)
        plt.tight_layout()
        plt.show()

        # Prepare CSV data
        for name, value in zip(names, values):
            csv_data.append({
                'Component': name,
                'Exergy_Destruction_W': value,
                'Exergy_Destruction_kW': value / 1e3,
                'Percentage': (value / total * 100) if total > 0 else 0
            })

        # Save CSV
        if save_csv and csv_data:
            df = pd.DataFrame(csv_data)
            csv_filename = self._generate_csv_filename('exergy_pie', csv_path)
            self._save_dataframe_to_csv(df, csv_filename, 'Exergy pie chart')

        return fig, ax

    # ================================================================ #
    #  7.  Energy Flow Summary Bar
    # ================================================================ #
    def plot_energy_summary(self, figsize=(10, 6), save_csv=False, csv_path=None):
        """
        Grouped bar chart: power in, power out, heat in, heat out.

        Parameters
        ----------
        figsize : tuple
            Figure size.
        save_csv : bool
            If True, save data to CSV.
        csv_path : str, optional
            Custom CSV file path.

        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        try:
            from .turbomachinery import Turbine, Pump, Compressor
            from .heat_exchangers import HeatExchanger
        except ImportError:
            # Fallback: try to identify components by attributes
            warnings.warn("Could not import component classes, using attribute detection")
            # Compressor was omitted here, so the isinstance() test below
            # raised NameError instead of degrading gracefully.
            Turbine = Pump = Compressor = HeatExchanger = None

        categories = {'Turbine Work': 0, 'Pump Work': 0,
                       'Heat Added': 0, 'Heat Rejected': 0}
        
        component_details = []

        for comp_id, comp in self.model.Component.items():
            comp_type = type(comp).__name__
            
            if Turbine and isinstance(comp, Turbine) and comp.work is not None:
                categories['Turbine Work'] += comp.work
                component_details.append({
                    'Component': comp_id,
                    'Type': 'Turbine',
                    'Value_W': comp.work,
                    'Value_kW': comp.work / 1e3
                })
            elif (((Pump and isinstance(comp, Pump))
                   or (Compressor and isinstance(comp, Compressor)))
                  and comp.work is not None):
                categories['Pump Work'] += comp.work
                component_details.append({
                    'Component': comp_id,
                    'Type': 'Pump',
                    'Value_W': comp.work,
                    'Value_kW': comp.work / 1e3
                })
            elif HeatExchanger and isinstance(comp, HeatExchanger):
                if (hasattr(comp, 'HeatAdded') and getattr(comp, 'Q', None)
                        is not None):
                    if comp.HeatAdded is True:
                        categories['Heat Added'] += comp.Q
                        component_details.append({
                            'Component': comp_id,
                            'Type': 'Heat_Added',
                            'Value_W': comp.Q,
                            'Value_kW': comp.Q / 1e3
                        })
                    elif comp.HeatAdded is False:
                        categories['Heat Rejected'] += comp.Q
                        component_details.append({
                            'Component': comp_id,
                            'Type': 'Heat_Rejected',
                            'Value_W': comp.Q,
                            'Value_kW': comp.Q / 1e3
                        })
            else:
                # Fallback: use attributes
                if getattr(comp, 'work', None) is not None:
                    if 'turbine' in comp_id.lower():
                        categories['Turbine Work'] += comp.work
                    elif 'pump' in comp_id.lower():
                        categories['Pump Work'] += comp.work

        fig, ax = plt.subplots(figsize=figsize)
        colors = ['green', 'red', 'orange', 'blue']
        bars = ax.bar(categories.keys(),
                      [v / 1e3 for v in categories.values()],
                      color=colors, edgecolor='black', linewidth=0.5)

        max_val = max(categories.values()) if categories.values() else 1
        
        for bar, v in zip(bars, categories.values()):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max_val/1e3*0.01,
                    f'{v/1e3:.2f} kW', ha='center', fontsize=11)

        ax.set_ylabel('Energy  [kW]', fontsize=13)
        ax.set_title('Energy Flow Summary', fontsize=15)
        plt.tight_layout()
        plt.show()

        # Save CSV
        if save_csv:
            # Create summary CSV
            summary_data = [{
                'Category': k,
                'Value_W': v,
                'Value_kW': v / 1e3
            } for k, v in categories.items()]
            
            df_summary = pd.DataFrame(summary_data)
            csv_filename = self._generate_csv_filename('energy_summary', csv_path)
            self._save_dataframe_to_csv(df_summary, csv_filename, 'Energy summary')
            
            # Also save detailed breakdown if available
            if component_details:
                df_details = pd.DataFrame(component_details)
                detail_filename = csv_filename.replace('.csv', '_detailed.csv')
                self._save_dataframe_to_csv(df_details, detail_filename, 'Energy summary (detailed)')

        return fig, ax

    # ================================================================ #
    #  8.  Export All State Points to CSV
    # ================================================================ #
    def export_all_points(self, csv_path=None):
        """
        Export all state points with all properties to CSV.

        Parameters
        ----------
        csv_path : str, optional
            Custom CSV file path. Auto-generates if None.

        Returns
        -------
        csv_path : str
            Path to saved CSV file.
        """
        data = []
        
        for name, pt in self.model.Point.items():
            row = {'Point': name}
            
            # Common properties
            props = ['T', 'P', 'H', 'S', 'D', 'Q', 'fluid']
            
            for prop in props:
                value = self._safe_get_property(pt, prop)
                row[prop] = value
            
            # Add converted units
            if row['T'] is not None:
                row['T_C'] = row['T'] - 273.15
            if row['P'] is not None:
                row['P_bar'] = row['P'] / 1e5
            if row['H'] is not None:
                row['H_kJ_per_kg'] = row['H'] / 1e3
            if row['S'] is not None:
                row['S_kJ_per_kg_K'] = row['S'] / 1e3
            
            data.append(row)
        
        if not data:
            warnings.warn("No state points to export")
            return None
        
        df = pd.DataFrame(data)
        csv_filename = self._generate_csv_filename('all_state_points', csv_path)
        
        if self._save_dataframe_to_csv(df, csv_filename, 'All state points'):
            return csv_filename
        return None

    # ================================================================ #
    #  9.  Export All Components to CSV
    # ================================================================ #
    def export_all_components(self, csv_path=None):
        """
        Export all component data to CSV.

        Parameters
        ----------
        csv_path : str, optional
            Custom CSV file path.

        Returns
        -------
        csv_path : str
            Path to saved CSV file.
        """
        data = []
        
        for comp_id, comp in self.model.Component.items():
            row = {
                'Component': comp_id,
                'Type': type(comp).__name__,
            }
            
            # Common attributes
            attrs = ['work', 'Q', 'Ex_D', 'eta', 'Solution_Status']
            
            for attr in attrs:
                if hasattr(comp, attr):
                    value = getattr(comp, attr)
                    row[attr] = value
                    
                    # Add converted units
                    if attr in ['work', 'Q', 'Ex_D'] and isinstance(value, (int, float)):
                        row[f'{attr}_kW'] = value / 1e3
            
            data.append(row)
        
        if not data:
            warnings.warn("No components to export")
            return None
        
        df = pd.DataFrame(data)
        csv_filename = self._generate_csv_filename('all_components', csv_path)
        
        if self._save_dataframe_to_csv(df, csv_filename, 'All components'):
            return csv_filename
        return None

    # ================================================================ #
    #  Internal:  Saturation Domes
    # ================================================================ #
    def _draw_saturation_dome(self, ax, fluid):
        """Draw saturation dome on T-s axes."""
        try:
            T_min = CP.PropsSI('TTRIPLE', fluid) + 1
            T_max = CP.PropsSI('TCRIT', fluid) - 1
            T_range = np.linspace(T_min, T_max, 300)

            S_liq = np.array([CP.PropsSI('S', 'T', T, 'Q', 0, fluid)
                              for T in T_range]) / 1e3
            S_vap = np.array([CP.PropsSI('S', 'T', T, 'Q', 1, fluid)
                              for T in T_range]) / 1e3
            T_plot = T_range - 273.15

            ax.plot(S_liq, T_plot, 'k-', linewidth=1, alpha=0.5, label='Saturation')
            ax.plot(S_vap, T_plot, 'k-', linewidth=1, alpha=0.5)
            ax.fill_betweenx(T_plot, S_liq, S_vap,
                             alpha=0.06, color='blue')
        except Exception as e:
            warnings.warn(f"Could not draw T-s saturation dome for {fluid}: {e}")

    def _draw_saturation_dome_Ph(self, ax, fluid):
        """Draw saturation dome on P-h axes."""
        try:
            T_min = CP.PropsSI('TTRIPLE', fluid) + 1
            T_max = CP.PropsSI('TCRIT', fluid) - 1
            T_range = np.linspace(T_min, T_max, 300)

            H_liq = np.array([CP.PropsSI('H', 'T', T, 'Q', 0, fluid)
                              for T in T_range]) / 1e3
            H_vap = np.array([CP.PropsSI('H', 'T', T, 'Q', 1, fluid)
                              for T in T_range]) / 1e3
            P_liq = np.array([CP.PropsSI('P', 'T', T, 'Q', 0, fluid)
                              for T in T_range]) / 1e5
            P_vap = np.array([CP.PropsSI('P', 'T', T, 'Q', 1, fluid)
                              for T in T_range]) / 1e5

            ax.semilogy(H_liq, P_liq, 'k-', linewidth=1, alpha=0.5, label='Saturation')
            ax.semilogy(H_vap, P_vap, 'k-', linewidth=1, alpha=0.5)
        except Exception as e:
            warnings.warn(f"Could not draw P-h saturation dome for {fluid}: {e}")

    def _draw_saturation_dome_hs(self, ax, fluid):
        """Draw saturation dome on h-s axes."""
        try:
            T_min = CP.PropsSI('TTRIPLE', fluid) + 1
            T_max = CP.PropsSI('TCRIT', fluid) - 1
            T_range = np.linspace(T_min, T_max, 300)

            S_liq = np.array([CP.PropsSI('S', 'T', T, 'Q', 0, fluid)
                              for T in T_range]) / 1e3
            S_vap = np.array([CP.PropsSI('S', 'T', T, 'Q', 1, fluid)
                              for T in T_range]) / 1e3
            H_liq = np.array([CP.PropsSI('H', 'T', T, 'Q', 0, fluid)
                              for T in T_range]) / 1e3
            H_vap = np.array([CP.PropsSI('H', 'T', T, 'Q', 1, fluid)
                              for T in T_range]) / 1e3

            ax.plot(S_liq, H_liq, 'k-', linewidth=1, alpha=0.5, label='Saturation')
            ax.plot(S_vap, H_vap, 'k-', linewidth=1, alpha=0.5)
        except Exception as e:
            warnings.warn(f"Could not draw h-s saturation dome for {fluid}: {e}")