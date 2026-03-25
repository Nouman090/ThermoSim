"""
plotting.py
-----------
All visualisation lives here — completely separated from computation.

Usage:
    from thermocycle.plotting import CyclePlotter
    plotter = CyclePlotter(Model)
    plotter.plot_Ts_diagram(['1','2','3','4','1'])
    plotter.plot_Ph_diagram(['1','2','3','4','1'])
    plotter.plot_hex_profile(Model.Component['HEX1'])
    plotter.plot_exergy_bar()
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import warnings


class CyclePlotter:
    """Handles every kind of plot for a ThermodynamicModel."""

    def __init__(self, model):
        """
        Parameters
        ----------
        model : ThermodynamicModel
            A solved model whose Point and Component dicts are populated.
        """
        self.model = model

    # ================================================================ #
    #  1.  T-s Diagram
    # ================================================================ #
    def plot_Ts_diagram(self, loop_points=None, fluid=None,
                        show_dome=True, show_labels=True,
                        title='T-s Diagram', figsize=(10, 7)):
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
        """
        fig, ax = plt.subplots(figsize=figsize)

        # ---- determine fluid for the dome ----------------------------
        if fluid is None and loop_points is not None:
            fluid = self.model.Point[loop_points[0]].fluid

        # ---- saturation dome -----------------------------------------
        if show_dome and fluid and fluid != 'Therminol66':
            self._draw_saturation_dome(ax, fluid)

        # ---- cycle line or scatter -----------------------------------
        if loop_points is not None:
            S_vals = []
            T_vals = []
            for name in loop_points:
                pt = self.model.Point[name]
                if pt.S is None or pt.T is None:
                    warnings.warn(f"Point '{name}' has no S or T — skipped.")
                    continue
                S_vals.append(pt.S / 1e3)           # kJ/(kg·K)
                T_vals.append(pt.T - 273.15)         # °C

            ax.plot(S_vals, T_vals, 'bo-',
                    linewidth=2, markersize=8, zorder=5)

            if show_labels:
                # label every point except the duplicate closer
                for i, name in enumerate(loop_points[:-1]):
                    ax.annotate(
                        name,
                        (S_vals[i], T_vals[i]),
                        textcoords='offset points',
                        xytext=(8, 8),
                        fontsize=12, fontweight='bold',
                        color='darkblue',
                    )
        else:
            # just scatter all points
            for name, pt in self.model.Point.items():
                if pt.S is not None and pt.T is not None:
                    ax.scatter(pt.S / 1e3, pt.T - 273.15,
                               s=60, zorder=5)
                    if show_labels:
                        ax.annotate(
                            name,
                            (pt.S / 1e3, pt.T - 273.15),
                            textcoords='offset points',
                            xytext=(5, 5), fontsize=10,
                        )

        ax.set_xlabel('Entropy  [kJ/(kg·K)]', fontsize=13)
        ax.set_ylabel('Temperature  [°C]', fontsize=13)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        return fig, ax

    # ================================================================ #
    #  2.  P-h Diagram
    # ================================================================ #
    def plot_Ph_diagram(self, loop_points=None, fluid=None,
                        show_dome=True, show_labels=True,
                        title='P-h Diagram', figsize=(10, 7)):
        """
        Plot a log-Pressure vs Enthalpy diagram.

        Parameters are the same as plot_Ts_diagram.
        """
        fig, ax = plt.subplots(figsize=figsize)

        if fluid is None and loop_points is not None:
            fluid = self.model.Point[loop_points[0]].fluid

        # dome
        if show_dome and fluid and fluid != 'Therminol66':
            self._draw_saturation_dome_Ph(ax, fluid)

        if loop_points is not None:
            H_vals = []
            P_vals = []
            for name in loop_points:
                pt = self.model.Point[name]
                if pt.H is None or pt.P is None:
                    continue
                H_vals.append(pt.H / 1e3)
                P_vals.append(pt.P / 1e5)

            ax.semilogy(H_vals, P_vals, 'ro-',
                        linewidth=2, markersize=8, zorder=5)

            if show_labels:
                for i, name in enumerate(loop_points[:-1]):
                    ax.annotate(
                        name,
                        (H_vals[i], P_vals[i]),
                        textcoords='offset points',
                        xytext=(8, 8),
                        fontsize=12, fontweight='bold',
                        color='darkred',
                    )
        else:
            for name, pt in self.model.Point.items():
                if pt.H is not None and pt.P is not None:
                    ax.scatter(pt.H / 1e3, pt.P / 1e5, s=60, zorder=5)
                    if show_labels:
                        ax.annotate(name, (pt.H/1e3, pt.P/1e5),
                                    textcoords='offset points',
                                    xytext=(5, 5), fontsize=10)

        ax.set_xlabel('Enthalpy  [kJ/kg]', fontsize=13)
        ax.set_ylabel('Pressure  [bar]', fontsize=13)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.3, which='both')
        plt.tight_layout()
        plt.show()
        return fig, ax

    # ================================================================ #
    #  3.  h-s Diagram  (Mollier)
    # ================================================================ #
    def plot_hs_diagram(self, loop_points=None, fluid=None,
                        show_dome=True, show_labels=True,
                        title='h-s Diagram', figsize=(10, 7)):
        """Enthalpy vs Entropy (Mollier) diagram."""
        fig, ax = plt.subplots(figsize=figsize)

        if fluid is None and loop_points is not None:
            fluid = self.model.Point[loop_points[0]].fluid

        if show_dome and fluid and fluid != 'Therminol66':
            self._draw_saturation_dome_hs(ax, fluid)

        if loop_points is not None:
            S_vals = [self.model.Point[n].S / 1e3 for n in loop_points
                      if self.model.Point[n].S is not None]
            H_vals = [self.model.Point[n].H / 1e3 for n in loop_points
                      if self.model.Point[n].H is not None]
            ax.plot(S_vals, H_vals, 'go-', linewidth=2, markersize=8, zorder=5)

            if show_labels:
                for i, name in enumerate(loop_points[:-1]):
                    ax.annotate(name, (S_vals[i], H_vals[i]),
                                textcoords='offset points',
                                xytext=(8, 8), fontsize=12,
                                fontweight='bold', color='darkgreen')

        ax.set_xlabel('Entropy  [kJ/(kg·K)]', fontsize=13)
        ax.set_ylabel('Enthalpy  [kJ/kg]', fontsize=13)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        return fig, ax

    # ================================================================ #
    #  4.  Heat-Exchanger Temperature Profile
    # ================================================================ #
    def plot_hex_profile(self, hex_id, div_N=200, figsize=(12, 5)):
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

        Returns
        -------
        min_dT : float
            The minimum temperature difference (pinch) in K.
        """
        hx = self.model.Component[hex_id]
        if not hx.Solution_Status:
            raise ValueError(f"{hex_id} has not been solved yet.")

        if hx.HEX_type == 'SimpleHEX':
            raise ValueError(
                "SimpleHEX has only two points — no internal profile."
            )

        Q_total = hx.Q
        q = Q_total / div_N
        dP_h = hx.Hot_In.P - hx.Hot_Out.P
        dP_c = hx.Cold_In.P - hx.Cold_Out.P

        Th = np.zeros(div_N + 1)
        Tc = np.zeros(div_N + 1)

        for n in range(div_N + 1):
            hh = hx.Hot_Out.H + q * n / hx.Hot_Mass_flowrate
            hc = hx.Cold_In.H + q * n / hx.Cold_Mass_flowrate
            Ph = hx.Hot_Out.P + (dP_h / div_N) * n
            Pc = hx.Cold_In.P - (dP_c / div_N) * n

            Th[n] = self.model.Prop(
                hx.Hot_In.fluid, StatePointName='_plot_h',
                H=hh, P=Ph
            ).T - 273.15
            Tc[n] = self.model.Prop(
                hx.Cold_In.fluid, StatePointName='_plot_c',
                H=hc, P=Pc
            ).T - 273.15

        dT = Th - Tc
        Q_axis = np.linspace(0, Q_total / 1e3, div_N + 1)

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
        ax2.axhline(y=min(dT), color='red', linestyle='--',
                     linewidth=1.5,
                     label=f'Min ΔT = {min(dT):.2f} K')
        idx_min = np.argmin(dT)
        ax2.plot(Q_axis[idx_min], dT[idx_min], 'rv', markersize=12)
        ax2.set_xlabel('Heat transferred  [kW]', fontsize=12)
        ax2.set_ylabel('ΔT  [K]', fontsize=12)
        ax2.set_title(f'{hex_id} — Pinch Diagram', fontsize=13)
        ax2.legend(fontsize=11)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()
        return min(dT)

    # ================================================================ #
    #  5.  Exergy Destruction Bar Chart
    # ================================================================ #
    def plot_exergy_bar(self, figsize=(12, 6), as_percentage=False):
        """
        Bar chart of exergy destruction for every component.

        Parameters
        ----------
        as_percentage : bool
            If True, show each component's share of total Ex_D.
        """
        names = []
        values = []
        for comp_id, comp in self.model.Component.items():
            if isinstance(comp.Ex_D, (int, float)):
                names.append(comp_id)
                values.append(comp.Ex_D)

        if not values:
            print("No numeric exergy destruction data found.")
            return

        values = np.array(values)
        total = values.sum()

        fig, ax = plt.subplots(figsize=figsize)
        colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(values)))

        if as_percentage and total > 0:
            pct = values / total * 100
            bars = ax.bar(names, pct, color=colors, edgecolor='black',
                          linewidth=0.5)
            ax.set_ylabel('Exergy Destruction  [%]', fontsize=13)

            # add value labels on bars
            for bar, p in zip(bars, pct):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.5,
                        f'{p:.1f}%', ha='center', fontsize=10)
        else:
            bars = ax.bar(names, values / 1e3, color=colors,
                          edgecolor='black', linewidth=0.5)
            ax.set_ylabel('Exergy Destruction  [kW]', fontsize=13)

            for bar, v in zip(bars, values / 1e3):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(values/1e3)*0.01,
                        f'{v:.2f}', ha='center', fontsize=10)

        ax.set_title('Component-wise Exergy Destruction', fontsize=15)
        plt.xticks(rotation=45, ha='right', fontsize=11)
        plt.tight_layout()
        plt.show()
        return fig, ax

    # ================================================================ #
    #  6.  Exergy Pie Chart
    # ================================================================ #
    def plot_exergy_pie(self, figsize=(8, 8)):
        """Pie chart showing share of exergy destruction."""
        names = []
        values = []
        for comp_id, comp in self.model.Component.items():
            if isinstance(comp.Ex_D, (int, float)) and comp.Ex_D > 0:
                names.append(comp_id)
                values.append(comp.Ex_D)

        if not values:
            print("No numeric exergy destruction data found.")
            return

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
        return fig, ax

    # ================================================================ #
    #  7.  Energy Flow Summary Bar
    # ================================================================ #
    def plot_energy_summary(self, figsize=(10, 6)):
        """Grouped bar chart: power in, power out, heat in, heat out."""
        from .turbomachinery import Turbine, Pump
        from .heat_exchangers import HeatExchanger

        categories = {'Turbine Work': 0, 'Pump Work': 0,
                       'Heat Added': 0, 'Heat Rejected': 0}

        for comp in self.model.Component.values():
            if isinstance(comp, Turbine):
                categories['Turbine Work'] += comp.work
            elif isinstance(comp, Pump):
                categories['Pump Work'] += comp.work
            elif isinstance(comp, HeatExchanger):
                if comp.HeatAdded is True:
                    categories['Heat Added'] += comp.Q
                elif comp.HeatAdded is False:
                    categories['Heat Rejected'] += comp.Q

        fig, ax = plt.subplots(figsize=figsize)
        colors = ['green', 'red', 'orange', 'blue']
        bars = ax.bar(categories.keys(),
                      [v / 1e3 for v in categories.values()],
                      color=colors, edgecolor='black', linewidth=0.5)

        for bar, v in zip(bars, categories.values()):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(categories.values())/1e3*0.01,
                    f'{v/1e3:.2f} kW', ha='center', fontsize=11)

        ax.set_ylabel('Energy  [kW]', fontsize=13)
        ax.set_title('Energy Flow Summary', fontsize=15)
        plt.tight_layout()
        plt.show()
        return fig, ax

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

            ax.plot(S_liq, T_plot, 'k-', linewidth=1, alpha=0.5)
            ax.plot(S_vap, T_plot, 'k-', linewidth=1, alpha=0.5)
            ax.fill_betweenx(T_plot, S_liq, S_vap,
                             alpha=0.06, color='blue')
        except Exception:
            pass

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

            ax.semilogy(H_liq, P_liq, 'k-', linewidth=1, alpha=0.5)
            ax.semilogy(H_vap, P_vap, 'k-', linewidth=1, alpha=0.5)
        except Exception:
            pass

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

            ax.plot(S_liq, H_liq, 'k-', linewidth=1, alpha=0.5)
            ax.plot(S_vap, H_vap, 'k-', linewidth=1, alpha=0.5)
        except Exception:
            pass