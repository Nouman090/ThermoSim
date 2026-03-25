"""
simple_components.py
--------------------
Source, Sink, Mixer, Splitter, Separator, Pipe, Expansion_valve
"""

import warnings
from .base_component import Component


# ================================================================== #
#   SOURCE
# ================================================================== #
class Source(Component):

    def __init__(self, Model, ID, Out_state, Calculate=False):
        self.Out_state = Out_state
        self.energy_supply = None
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.Out = self.Model.Point[self.Out_state]
        if self.Out.Mass_flowrate is None or self.Out.H is None:
            raise ValueError(
                f"Mass-flow rate or enthalpy of {self.ID} outlet is missing."
            )
        self.energy_supply = self.Out.H * self.Out.Mass_flowrate
        self.Ex_D = 0
        self.Solution_Status = True

    def __str__(self):
        return (
            f"{self.ID} (Source):\n"
            f"  Energy supply  : {self.energy_supply} W\n"
            f"  Solved         : {self.Solution_Status}"
        )


# ================================================================== #
#   SINK
# ================================================================== #
class Sink(Component):

    def __init__(self, Model, ID, In_state, Calculate=False):
        self.In_state = In_state
        self.energy_supply = None
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In = self.Model.Point[self.In_state]
        if self.In.Mass_flowrate is None or self.In.H is None:
            raise ValueError(
                f"Mass-flow rate or enthalpy of {self.ID} inlet is missing."
            )
        self.energy_supply = self.In.H * self.In.Mass_flowrate
        self.Ex_D = 0
        self.Solution_Status = True

    def __str__(self):
        return (
            f"{self.ID} (Sink):\n"
            f"  Energy supply  : {self.energy_supply} W\n"
            f"  Solved         : {self.Solution_Status}"
        )


# ================================================================== #
#   MIXER
# ================================================================== #
class Mixer(Component):

    def __init__(self, Model, ID, In_states, Out_state, Calculate=False):
        self.In_states = In_states        # list of state-point names
        self.Out_state = Out_state
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.Out = self.Model.Point[self.Out_state]

        mass_in   = 0.0
        energy_in = 0.0
        for name in self.In_states:
            pt = self.Model.Point[name]
            if pt.Mass_flowrate is None or pt.H is None:
                raise ValueError(
                    f"Mass-flow rate or enthalpy for '{name}' "
                    f"in Mixer {self.ID} is missing."
                )
            mass_in   += pt.Mass_flowrate
            energy_in += pt.Mass_flowrate * pt.H

        h_out = energy_in / mass_in

        # keep existing outlet pressure, else take first inlet's
        if self.Out.P is None:
            self.Out.P = self.Model.Point[self.In_states[0]].P

        # if outlet enthalpy was already set, check consistency
        if self.Out.H is not None:
            if abs(self.Out.H - h_out) > 1.0:      # 1 J/kg tolerance
                raise ValueError(
                    f"Mixer {self.ID}: given outlet H = {self.Out.H:.2f} "
                    f"but energy balance gives {h_out:.2f}"
                )

        self.Out = self.Model.Prop(
            self.Out.fluid,
            StatePointName=self.Out.StatePointName,
            P=self.Out.P, H=h_out,
            Mass_flowrate=mass_in
        )

        # exergy destruction
        try:
            ex_in  = sum(self.Model.Point[s].Ex for s in self.In_states)
            ex_out = self.Out.Ex
            self.Ex_D = ex_in - ex_out
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True
        self.Model.Point[self.Out_state] = self.Out

    def __str__(self):
        return (
            f"{self.ID} (Mixer):\n"
            f"  Total mass flow : {self.Out.Mass_flowrate:.4f} kg/s\n"
            f"  Outlet H        : {self.Out.H:.2f} J/kg\n"
            f"  Exergy destr.   : {self.Ex_D}\n"
            f"  Solved          : {self.Solution_Status}"
        )


# ================================================================== #
#   SPLITTER
# ================================================================== #
class Splitter(Component):

    def __init__(self, Model, ID, In_state, Out_states,
                 split_fractions, Calculate=False):
        self.In_state = In_state
        self.Out_states = Out_states
        self.split_fractions = list(split_fractions)
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In = self.Model.Point[self.In_state]

        if self.In.Mass_flowrate is None or self.In.H is None:
            raise ValueError(
                f"Inlet data for Splitter {self.ID} is incomplete."
            )

        # normalise fractions if they don't sum to 1
        total = sum(self.split_fractions)
        if abs(total - 1.0) > 1e-5:
            warnings.warn(
                f"Split fractions for {self.ID} sum to {total}. "
                f"Normalising."
            )
            self.split_fractions = [f / total for f in self.split_fractions]

        for i, name in enumerate(self.Out_states):
            m_out = self.In.Mass_flowrate * self.split_fractions[i]
            pt = self.Model.Prop(
                self.In.fluid,
                StatePointName=name,
                P=self.In.P, H=self.In.H,
                Mass_flowrate=m_out
            )
            self.Model.Point[name] = pt

        self.Ex_D = 0
        self.Solution_Status = True

    def __str__(self):
        return (
            f"{self.ID} (Splitter):\n"
            f"  Inlet mass flow  : {self.In.Mass_flowrate:.4f} kg/s\n"
            f"  Split fractions  : {self.split_fractions}\n"
            f"  Exergy destr.    : {self.Ex_D}\n"
            f"  Solved           : {self.Solution_Status}"
        )


# ================================================================== #
#   SEPARATOR  (flash drum)
# ================================================================== #
class Separator(Component):

    def __init__(self, Model, ID, In_state,
                 Out_vap_state, Out_liq_state, Calculate=False):
        self.In_state = In_state
        self.Out_vap_state = Out_vap_state
        self.Out_liq_state = Out_liq_state
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In      = self.Model.Point[self.In_state]
        self.Out_vap = self.Model.Point[self.Out_vap_state]
        self.Out_liq = self.Model.Point[self.Out_liq_state]

        quality = self.In.Q
        if quality is None or not (0 <= quality <= 1):
            warnings.warn(
                f"[{self.ID}] Inlet '{self.In_state}' is not two-phase "
                f"(Q = {quality})."
            )

        P = self.In.P
        m_vap = self.In.Mass_flowrate * quality
        m_liq = self.In.Mass_flowrate * (1.0 - quality)

        # saturated vapour outlet
        self.Out_vap = self.Model.Prop(
            self.In.fluid, StatePointName=self.Out_vap_state,
            P=P, Q=1.0, Mass_flowrate=m_vap
        )
        # saturated liquid outlet
        self.Out_liq = self.Model.Prop(
            self.In.fluid, StatePointName=self.Out_liq_state,
            P=P, Q=0.0, Mass_flowrate=m_liq
        )

        try:
            self.Ex_D = self.In.Ex - self.Out_vap.Ex - self.Out_liq.Ex
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True
        self.Model.Point[self.Out_vap_state] = self.Out_vap
        self.Model.Point[self.Out_liq_state] = self.Out_liq

    def __str__(self):
        return (
            f"{self.ID} (Separator):\n"
            f"  Inlet flow       : {self.In.Mass_flowrate:.2f} kg/s\n"
            f"  Vapour outlet    : {self.Out_vap.Mass_flowrate:.2f} kg/s\n"
            f"  Liquid outlet    : {self.Out_liq.Mass_flowrate:.2f} kg/s\n"
            f"  Exergy destr.    : {self.Ex_D}\n"
            f"  Solved           : {self.Solution_Status}"
        )


# ================================================================== #
#   PIPE
# ================================================================== #
class Pipe(Component):

    def __init__(self, Model, ID, In_state, Out_state,
                 Pressure_drop=0, Temperature_drop=0, Calculate=False):
        self.In_state = In_state
        self.Out_state = Out_state
        self.Pressure_drop = Pressure_drop
        self.Temperature_drop = Temperature_drop
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]
        _, _ = self._resolve_mass_flowrate(self.In, self.Out)

        if self.In.P is None:
            self.In.P = self.Out.P + self.Pressure_drop
            self.In.T = self.Out.T + self.Temperature_drop
        elif self.Out.P is None:
            self.Out.P = self.In.P - self.Pressure_drop
            self.Out.T = self.In.T - self.Temperature_drop
        else:
            self.Pressure_drop    = self.In.P - self.Out.P
            self.Temperature_drop = self.In.T - self.Out.T

        try:
            self.Ex_D = self.In.Ex - self.Out.Ex
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True
        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        return (
            f"{self.ID} (Pipe):\n"
            f"  ΔP  : {self.Pressure_drop/1e5:.2f} bar\n"
            f"  ΔT  : {self.Temperature_drop:.2f} K\n"
            f"  Ex_D: {self.Ex_D}\n"
            f"  Solved: {self.Solution_Status}"
        )


# ================================================================== #
#   EXPANSION VALVE
# ================================================================== #
class Expansion_valve(Component):

    def __init__(self, Model, ID, In_state, Out_state, Calculate=False):
        self.In_state = In_state
        self.Out_state = Out_state
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]
        _, _ = self._resolve_mass_flowrate(self.In, self.Out)

        # isenthalpic expansion
        if self.In.H is None and self.Out.H is not None:
            self.In = self.Model.Prop(
                self.In.fluid,
                StatePointName=self.In.StatePointName,
                P=self.In.P, H=self.Out.H
            )
        elif self.In.H is not None and self.Out.H is None:
            self.Out = self.Model.Prop(
                self.Out.fluid,
                StatePointName=self.Out.StatePointName,
                P=self.Out.P, H=self.In.H
            )
        elif self.In.H is not None and self.Out.H is not None:
            if abs(self.In.H - self.Out.H) > 1.0:
                raise ValueError(
                    f"Enthalpy mismatch across {self.ID}: "
                    f"In.H={self.In.H:.2f}, Out.H={self.Out.H:.2f}"
                )
        else:
            raise ValueError(
                f"At least one side of {self.ID} needs a known enthalpy."
            )

        try:
            self.Ex_D = self.In.Ex - self.Out.Ex
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True
        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        return (
            f"{self.ID} (Expansion valve):\n"
            f"  P_in  : {self.In.P/1e5:.2f} bar\n"
            f"  P_out : {self.Out.P/1e5:.2f} bar\n"
            f"  Ex_D  : {self.Ex_D}\n"
            f"  Solved: {self.Solution_Status}"
        )