"""
simple_components.py
--------------------
Source, Sink, Mixer, Splitter, Separator, Pipe, Expansion_valve
"""

import warnings
from .base_component import Component


# ================================================================== #
#   VALIDATION DECORATOR
# ================================================================== #
def validate_inputs(func):
    """Decorator to check component inputs before calculation"""
    def wrapper(self):
        # Check if component is part of a model
        if not hasattr(self, 'Model'):
            raise RuntimeError(f"{self.ID}: Not attached to a model")
        
        # Check state points exist
        for attr in dir(self):
            if attr.endswith('_state') or attr.endswith('_states'):
                states = getattr(self, attr)
                if isinstance(states, list):
                    for s in states:
                        if s not in self.Model.Point:
                            raise KeyError(f"{self.ID}: State '{s}' not found in model")
                elif states and states not in self.Model.Point:
                    raise KeyError(f"{self.ID}: State '{states}' not found in model")
        
        return func(self)
    return wrapper


# ================================================================== #
#   SOURCE
# ================================================================== #
class Source(Component):
    """Energy source component - adds fluid to the system"""

    def __init__(self, Model, ID, Out_state, Calculate=False):
        self.Out_state = Out_state
        self.energy_supply = None
        super().__init__(Model, ID, Calculate)

    @validate_inputs
    def Cal(self):
        self.Out = self.Model.Point[self.Out_state]
        
        # Validation
        if self.Out.Mass_flowrate is None:
            raise ValueError(f"{self.ID}: Outlet mass flow rate is missing")
        if self.Out.H is None:
            raise ValueError(f"{self.ID}: Outlet enthalpy is missing")
        if self.Out.Mass_flowrate <= 0:
            raise ValueError(f"{self.ID}: Mass flow rate must be positive")
        
        self.energy_supply = self.Out.H * self.Out.Mass_flowrate
        
        # Exergy destruction (negative for source)
        if self.Out.Ex is not None:
            self.Ex_D = self.Out.Ex
        else:
            self.Ex_D = None
            
        self.Solution_Status = True

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Source): not yet solved."
        return (
            f"{self.ID} (Source):\n"
            f"  Energy supply  : {self.energy_supply:.2f} W\n"
            f"  Mass flow      : {self.Out.Mass_flowrate:.4f} kg/s\n"
            f"  Exergy         : {self.Ex_D}\n"
            f"  Solved         : {self.Solution_Status}"
        )


# ================================================================== #
#   SINK
# ================================================================== #
class Sink(Component):
    """Energy sink component - removes fluid from the system"""

    def __init__(self, Model, ID, In_state, Calculate=False):
        self.In_state = In_state
        self.energy_supply = None
        super().__init__(Model, ID, Calculate)

    @validate_inputs
    def Cal(self):
        self.In = self.Model.Point[self.In_state]
        
        # Validation
        if self.In.Mass_flowrate is None:
            raise ValueError(f"{self.ID}: Inlet mass flow rate is missing")
        if self.In.H is None:
            raise ValueError(f"{self.ID}: Inlet enthalpy is missing")
        if self.In.Mass_flowrate <= 0:
            raise ValueError(f"{self.ID}: Mass flow rate must be positive")
        
        self.energy_supply = self.In.H * self.In.Mass_flowrate
        
        # Exergy destruction (negative for sink)
        if self.In.Ex is not None:
            self.Ex_D = -self.In.Ex
        else:
            self.Ex_D = None
            
        self.Solution_Status = True

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Sink): not yet solved."
        return (
            f"{self.ID} (Sink):\n"
            f"  Energy supply  : {self.energy_supply:.2f} W\n"
            f"  Mass flow      : {self.In.Mass_flowrate:.4f} kg/s\n"
            f"  Exergy         : {self.Ex_D}\n"
            f"  Solved         : {self.Solution_Status}"
        )


# ================================================================== #
#   MIXER
# ================================================================== #
class Mixer(Component):
    """Mixes multiple inlet streams into one outlet stream"""
    
    TOLERANCE_ENTHALPY = 1.0e-3      # J/kg
    PRESSURE_WARNING_THRESHOLD = 0.05  # 5% pressure difference triggers warning

    def __init__(self, Model, ID, In_states, Out_state, 
                 pressure_loss_fraction=0.0, Calculate=False):
        """
        Args:
            pressure_loss_fraction: Fractional pressure loss (e.g., 0.02 = 2% loss)
        """
        if not In_states or len(In_states) < 2:
            raise ValueError(f"Mixer {ID}: Requires at least 2 inlet states")
            
        self.In_states = In_states
        self.Out_state = Out_state
        self.pressure_loss_fraction = pressure_loss_fraction
        super().__init__(Model, ID, Calculate)

    @validate_inputs
    def Cal(self):
        self.Out = self.Model.Point[self.Out_state]
        inlets = {n: self.Model.Point[n] for n in self.In_states}

        # ── Classify ports ─────────────────────────────────────────────────
        # A port is "flow-known" when both Mass_flowrate and H are set,
        # because those two are the minimum needed to contribute to the
        # mass and energy balances.
        def flow_known(pt):
            return pt.Mass_flowrate is not None and pt.H is not None

        known_inlet_names   = [n for n, pt in inlets.items() if     flow_known(pt)]
        unknown_inlet_names = [n for n, pt in inlets.items() if not flow_known(pt)]
        outlet_known        = flow_known(self.Out)

        n_unknown_inlets = len(unknown_inlet_names)

        # ── Guard: only one unknown port is solvable ───────────────────────
        total_unknowns = n_unknown_inlets + (0 if outlet_known else 1)
        if total_unknowns > 1:
            raise ValueError(
                f"Mixer {self.ID}: {total_unknowns} ports are missing Mass_flowrate "
                f"or H. At most 1 port may be unknown (the others must be fully set)."
            )
        if total_unknowns == 0:
            # All ports already defined — just validate consistency and rebuild
            pass  # falls through to the forward solve below

        # ── Case A: one inlet is unknown → back-calculate it ──────────────
        if n_unknown_inlets == 1 and outlet_known:
            unk_name = unknown_inlet_names[0]
            unk_pt   = inlets[unk_name]

            # Validate outlet
            if self.Out.Mass_flowrate is None or self.Out.Mass_flowrate <= 0:
                raise ValueError(
                    f"Mixer {self.ID}: Outlet mass flow must be known and positive "
                    f"to back-calculate inlet '{unk_name}'."
                )
            if self.Out.P is None:
                raise ValueError(
                    f"Mixer {self.ID}: Outlet pressure must be known to "
                    f"back-calculate inlet '{unk_name}'."
                )

            # Sum contributions of the known inlets
            mass_known   = sum(inlets[n].Mass_flowrate for n in known_inlet_names)
            energy_known = sum(inlets[n].Mass_flowrate * inlets[n].H
                               for n in known_inlet_names)

            m_unk = self.Out.Mass_flowrate - mass_known
            if m_unk <= 0:
                raise ValueError(
                    f"Mixer {self.ID}: Back-calculated mass flow for inlet "
                    f"'{unk_name}' is non-positive ({m_unk:.6f} kg/s). "
                    f"Check whether the outlet mass flow exceeds the sum of "
                    f"the known inlets."
                )

            h_unk = (self.Out.Mass_flowrate * self.Out.H - energy_known) / m_unk

            # Write the derived state back into the model
            unk_pt.Mass_flowrate = m_unk
            unk_pt.H             = h_unk
            if unk_pt.P is None:
                unk_pt.P = self.Out.P   # best available estimate
            self.Model.Point[unk_name] = unk_pt
            # Refresh local dict so the forward pass sees the updated values
            inlets[unk_name] = unk_pt

        # ── Forward solve: all inlets now known, compute / validate outlet ─
        mass_in         = 0.0
        energy_in       = 0.0
        inlet_pressures = []

        for name, pt in inlets.items():
            if pt.Mass_flowrate is None or pt.H is None:
                raise ValueError(
                    f"Mixer {self.ID}: Mass flow or enthalpy still missing for "
                    f"inlet '{name}' after back-calculation attempt."
                )
            if pt.Mass_flowrate < 0:
                raise ValueError(
                    f"Mixer {self.ID}: Negative mass flow for inlet '{name}'."
                )
            if pt.P is None:
                raise ValueError(
                    f"Mixer {self.ID}: Pressure missing for inlet '{name}'."
                )
            mass_in   += pt.Mass_flowrate
            energy_in += pt.Mass_flowrate * pt.H
            inlet_pressures.append(pt.P)

        if mass_in <= 0:
            raise ValueError(
                f"Mixer {self.ID}: Total inlet mass flow must be positive."
            )

        h_out = energy_in / mass_in

        # Pressure handling
        P_min = min(inlet_pressures)
        P_max = max(inlet_pressures)
        if (P_max - P_min) / P_max > self.PRESSURE_WARNING_THRESHOLD:
            warnings.warn(
                f"Mixer {self.ID}: Inlet pressures vary significantly "
                f"({P_min/1e5:.2f} – {P_max/1e5:.2f} bar). Using minimum."
            )
        P_out = P_min * (1.0 - self.pressure_loss_fraction)

        if self.Out.P is None:
            self.Out.P = P_out
        elif abs(self.Out.P - P_out) / P_out > 0.001:
            warnings.warn(
                f"Mixer {self.ID}: Specified outlet pressure ({self.Out.P/1e5:.2f} bar) "
                f"differs from calculated ({P_out/1e5:.2f} bar)."
            )

        # Validate outlet enthalpy if already set
        if self.Out.H is not None:
            if abs(self.Out.H - h_out) > self.TOLERANCE_ENTHALPY:
                raise ValueError(
                    f"Mixer {self.ID}: Outlet enthalpy mismatch. "
                    f"Given: {self.Out.H:.2f} J/kg, Calculated: {h_out:.2f} J/kg"
                )

        # Create / update outlet state point
        self.Out = self.Model.Prop(
            self.Out.fluid,
            StatePointName=self.Out.StatePointName,
            P=self.Out.P, H=h_out,
            Mass_flowrate=mass_in
        )

        # Calculate exergy destruction
        try:
            ex_in = sum(self.Model.Point[s].Ex for s in self.In_states)
            self.Ex_D = ex_in - self.Out.Ex
        except (AttributeError, TypeError):
            self.Ex_D = None

        self.Solution_Status = True
        self.Model.Point[self.Out_state] = self.Out

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Mixer): not yet solved."
        inlet_flows = ", ".join([
            ("?" if self.Model.Point[s].Mass_flowrate is None
             else f"{self.Model.Point[s].Mass_flowrate:.4f}")
            for s in self.In_states
        ])
        return (
            f"{self.ID} (Mixer):\n"
            f"  Inlet flows    : [{inlet_flows}] kg/s\n"
            f"  Total outlet   : {self.Out.Mass_flowrate:.4f} kg/s\n"
            f"  Outlet H       : {self.Out.H:.2f} J/kg\n"
            f"  Outlet P       : {self.Out.P/1e5:.2f} bar\n"
            f"  Exergy destr.  : {self.Ex_D}\n"
            f"  Solved         : {self.Solution_Status}"
        )


# ================================================================== #
#   SPLITTER
# ================================================================== #
class Splitter(Component):
    """Splits one inlet stream into multiple outlet streams"""
    
    TOLERANCE_FRACTION_SUM = 1e-5

    def __init__(self, Model, ID, In_state, Out_states,
                 split_fractions, Calculate=False):
        """
        Args:
            split_fractions: List of fractions for each outlet (must sum to ~1.0)
        """
        if not Out_states or len(Out_states) < 2:
            raise ValueError(f"Splitter {ID}: Requires at least 2 outlet states")
            
        if len(Out_states) != len(split_fractions):
            raise ValueError(
                f"Splitter {ID}: Number of outlets ({len(Out_states)}) must match "
                f"number of split fractions ({len(split_fractions)})"
            )
        
        self.In_state = In_state
        self.Out_states = Out_states
        self.split_fractions = list(split_fractions)
        super().__init__(Model, ID, Calculate)

    @validate_inputs
    def Cal(self):
        self.In = self.Model.Point[self.In_state]

        # ── Validate and normalise fractions first ─────────────────────────
        for i, frac in enumerate(self.split_fractions):
            if frac < 0:
                raise ValueError(
                    f"Splitter {self.ID}: Negative split fraction at outlet {i}"
                )
            if frac > 1:
                raise ValueError(
                    f"Splitter {self.ID}: Split fraction > 1 at outlet {i}"
                )
        total = sum(self.split_fractions)
        if abs(total - 1.0) > self.TOLERANCE_FRACTION_SUM:
            warnings.warn(
                f"Splitter {self.ID}: Split fractions sum to {total:.6f}, "
                f"normalizing to 1.0"
            )
            self.split_fractions = [f / total for f in self.split_fractions]

        # ── Step 1 – Resolve thermodynamic state (P and H) ────────────────
        # A splitter is isobaric and isenthalpic: every port carries the
        # same P and H.  We accept that state from *any* port — inlet or
        # any outlet — so the user can define conditions on whichever port
        # is convenient.
        all_port_names = [self.In_state] + list(self.Out_states)
        all_points     = [self.Model.Point[n] for n in all_port_names]

        ref_P = next((pt.P for pt in all_points if pt.P is not None), None)
        ref_H = next((pt.H for pt in all_points if pt.H is not None), None)

        if ref_P is None:
            raise ValueError(
                f"Splitter {self.ID}: Pressure (P) is not set on any port. "
                f"Set P on at least one inlet or outlet."
            )
        if ref_H is None:
            raise ValueError(
                f"Splitter {self.ID}: Enthalpy (H) is not set on any port. "
                f"Set H on at least one inlet or outlet."
            )

        # Propagate P and H to every port that is still missing them
        for pt in all_points:
            if pt.P is None:
                pt.P = ref_P
            if pt.H is None:
                pt.H = ref_H

        # Re-fetch inlet in case it was just updated
        self.In = self.Model.Point[self.In_state]

        # ── Step 2 – Resolve mass flow rate ───────────────────────────────
        # The INLET is authoritative whenever it is known:
        #
        #   m_outlet_i = m_in * fraction_i
        #
        # Only when the inlet flow is unknown do we infer it backwards from
        # an outlet.  Up to v3.2.1 every port was polled and the results were
        # averaged, which silently corrupted any model that re-ran a splitter
        # with new fractions (e.g. the regenerative Rankine cycle in main.py,
        # where the extraction fraction y is refined after a first pass): the
        # outlets still held flows computed from the OLD fractions, so the
        # mean of the stale outlets and the correct inlet was wrong, and mass
        # conservation broke downstream.
        if self.In.Mass_flowrate is not None:
            m_in = self.In.Mass_flowrate
        else:
            # Infer from whichever outlets are known; they must agree.
            back_calc = []
            for i, name in enumerate(self.Out_states):
                pt = self.Model.Point[name]
                if pt.Mass_flowrate is not None:
                    frac = self.split_fractions[i]
                    if frac <= 0:
                        raise ValueError(
                            f"Splitter {self.ID}: Cannot infer m_in from outlet "
                            f"'{name}' because its split fraction is zero."
                        )
                    back_calc.append((name, pt.Mass_flowrate / frac))

            if not back_calc:
                # No mass flow information at all — thermodynamic state is
                # fully propagated above; we just cannot set flow rates yet.
                warnings.warn(
                    f"Splitter {self.ID}: No mass flow rate is set on any port. "
                    f"Thermodynamic properties (P, H) have been propagated to all "
                    f"ports, but mass flow rates remain unresolved."
                )
                self.Ex_D = 0
                self.Solution_Status = True   # thermodynamic state is solved
                return

            values = [v for _, v in back_calc]
            if max(values) - min(values) > 1e-6 * max(values):
                detail = ", ".join(f"{src}→{v:.6f}" for src, v in back_calc)
                warnings.warn(
                    f"Splitter {self.ID}: Outlet mass flows imply different "
                    f"inlet flows ({detail}). Using the mean; check the split "
                    f"fractions."
                )
            m_in = sum(values) / len(values)

        if m_in <= 0:
            raise ValueError(
                f"Splitter {self.ID}: Resolved inlet mass flow must be positive "
                f"(got {m_in:.6f} kg/s)."
            )

        # ── Step 3 – Write fully-resolved state points back to the model ──
        # Inlet
        self.Model.Point[self.In_state] = self.Model.Prop(
            self.In.fluid,
            StatePointName=self.In_state,
            P=ref_P, H=ref_H,
            Mass_flowrate=m_in
        )
        self.In = self.Model.Point[self.In_state]

        # Outlets
        for i, name in enumerate(self.Out_states):
            m_out = m_in * self.split_fractions[i]
            self.Model.Point[name] = self.Model.Prop(
                self.In.fluid,
                StatePointName=name,
                P=ref_P, H=ref_H,
                Mass_flowrate=m_out
            )

        # No exergy destruction in ideal splitter
        self.Ex_D = 0
        self.Solution_Status = True

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Splitter): not yet solved."
        outlet_flows = ", ".join([
            f"{self.Model.Point[s].Mass_flowrate}" 
            for s in self.Out_states
        ])
        return (
            f"{self.ID} (Splitter):\n"
            f"  Inlet mass flow  : {self.In.Mass_flowrate} kg/s\n"
            f"  Split fractions  : {[f'{x:.4f}' for x in self.split_fractions]}\n"
            f"  Outlet flows     : [{outlet_flows}] kg/s\n"
            f"  Exergy destr.    : {self.Ex_D}\n"
            f"  Solved           : {self.Solution_Status}"
        )


# ================================================================== #
#   SEPARATOR  (Flash Drum)
# ================================================================== #
class Separator(Component):
    """Separates two-phase mixture into vapor and liquid streams"""

    def __init__(self, Model, ID, In_state,
                 Out_vap_state, Out_liq_state, Calculate=False):
        self.In_state = In_state
        self.Out_vap_state = Out_vap_state
        self.Out_liq_state = Out_liq_state
        super().__init__(Model, ID, Calculate)

    @validate_inputs
    def Cal(self):
        self.In      = self.Model.Point[self.In_state]
        self.Out_vap = self.Model.Point[self.Out_vap_state]
        self.Out_liq = self.Model.Point[self.Out_liq_state]

        # ── Tier 1: Resolve pressure ───────────────────────────────────────
        # A separator (flash drum) is isobaric: all three ports share the
        # same saturation pressure.  Accept P from whichever port has it.
        P = self.In.P or self.Out_vap.P or self.Out_liq.P
        if P is None:
            raise ValueError(
                f"Separator {self.ID}: Pressure is not set on any port. "
                f"Set P on the inlet, the vapour outlet, or the liquid outlet."
            )
        # Propagate to ports that are still missing P
        for pt in (self.In, self.Out_vap, self.Out_liq):
            if pt.P is None:
                pt.P = P

        # ── Tier 2: Propagate saturated outlet states ──────────────────────
        # The saturated vapour (Q=1) and liquid (Q=0) enthalpies at pressure
        # P are uniquely determined by P alone — no mass flow needed.
        # We pre-populate the outlet points so downstream components can
        # read their thermodynamic properties immediately.
        vap_probe = self.Model.Prop(
            self.In.fluid, StatePointName="_sep_vap_probe", P=P, Q=1.0
        )
        liq_probe = self.Model.Prop(
            self.In.fluid, StatePointName="_sep_liq_probe", P=P, Q=0.0
        )
        h_vap_sat = vap_probe.H
        h_liq_sat = liq_probe.H

        # Write saturation enthalpies into outlet points if not yet set
        if self.Out_vap.H is None:
            self.Out_vap.H = h_vap_sat
        if self.Out_liq.H is None:
            self.Out_liq.H = h_liq_sat

        # ── Tier 3: Mass-flow split (requires inlet quality) ──────────────
        # Quality Q is computable only when the inlet enthalpy H is known.
        # If it isn't, we leave mass flows unresolved but mark the
        # thermodynamic state as solved.
        # Q is numeric in [0, 1] inside the dome and None outside it.
        # (Up to v3.2.1 Prop._classify_phase replaced Q with a descriptive
        # string such as "Superheated (46.97 K)" for single-phase states, so
        # the comparisons below raised
        # `TypeError: '<' not supported between 'float' and 'str'`.)
        quality = self.In.Q

        if quality is None:
            if self.In.H is None:
                warnings.warn(
                    f"Separator {self.ID}: Inlet quality cannot be determined "
                    f"because the inlet enthalpy is not set. Outlet P and "
                    f"saturated H values have been propagated; the mass-flow "
                    f"split is deferred."
                )
                self.Solution_Status = True   # thermodyn. state is known
                return
            raise ValueError(
                f"Separator {self.ID}: the inlet is single-phase "
                f"({self.In.phase}) at {self.In.T - 273.15:.2f} °C, "
                f"{self.In.P / 1e5:.2f} bar. A flash drum needs a two-phase "
                f"inlet (0 < Q < 1)."
            )

        if not (0.0 < quality < 1.0):
            raise ValueError(
                f"Separator {self.ID}: Inlet quality = {quality:.4f}. "
                f"Must be strictly between 0 and 1 (two-phase region); "
                f"Q=0 and Q=1 are saturated states with nothing to separate."
            )

        # Mass-flow split requires the inlet mass flow rate
        m_total = self.In.Mass_flowrate
        if m_total is None:
            warnings.warn(
                f"Separator {self.ID}: Inlet mass flow rate is not set. "
                f"Outlet enthalpies are propagated; mass-flow split deferred."
            )
            self.Solution_Status = True
            return

        if m_total <= 0:
            raise ValueError(
                f"Separator {self.ID}: Inlet mass flow must be positive "
                f"(got {m_total:.6f} kg/s)."
            )

        m_vap = m_total * quality
        m_liq = m_total * (1.0 - quality)

        # Fully-specified outlet states
        self.Out_vap = self.Model.Prop(
            self.In.fluid,
            StatePointName=self.Out_vap_state,
            P=P, Q=1.0,
            Mass_flowrate=m_vap
        )
        self.Out_liq = self.Model.Prop(
            self.In.fluid,
            StatePointName=self.Out_liq_state,
            P=P, Q=0.0,
            Mass_flowrate=m_liq
        )

        # Energy balance sanity check
        h_in = self.In.H
        h_out_check = (m_vap * self.Out_vap.H + m_liq * self.Out_liq.H) / m_total
        if abs(h_in - h_out_check) > 1.0:   # 1 J/kg tolerance
            warnings.warn(
                f"Separator {self.ID}: Energy balance deviation = "
                f"{abs(h_in - h_out_check):.2f} J/kg."
            )

        # Exergy destruction
        try:
            self.Ex_D = self.In.Ex - (self.Out_vap.Ex + self.Out_liq.Ex)
        except (AttributeError, TypeError):
            self.Ex_D = None

        self.Solution_Status = True
        self.Model.Point[self.Out_vap_state] = self.Out_vap
        self.Model.Point[self.Out_liq_state] = self.Out_liq

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Separator): not yet solved."
        return (
            f"{self.ID} (Separator):\n"
            f"  Inlet flow       : {self.In.Mass_flowrate:.4f} kg/s (Q={self.In.Q:.4f})\n"
            f"  Vapor outlet     : {self.Out_vap.Mass_flowrate:.4f} kg/s\n"
            f"  Liquid outlet    : {self.Out_liq.Mass_flowrate:.4f} kg/s\n"
            f"  Pressure         : {self.In.P/1e5:.2f} bar\n"
            f"  Exergy destr.    : {self.Ex_D}\n"
            f"  Solved           : {self.Solution_Status}"
        )


# ================================================================== #
#   PIPE
# ================================================================== #
class Pipe(Component):
    """Pipe with pressure and temperature drops"""
    
    TOLERANCE_PRESSURE = 1e3      # Pa
    TOLERANCE_TEMPERATURE = 0.1   # K

    def __init__(self, Model, ID, In_state, Out_state,
                 Pressure_drop=None, Temperature_drop=None, Calculate=False):
        """
        Args:
            Pressure_drop: Pressure drop in Pa (positive = pressure decreases).
                None means "not specified" - derive it from the two end
                pressures.  0 means "explicitly lossless".
            Temperature_drop: Temperature drop in K (positive = temperature
                decreases).  None / 0 as above.
        """
        self.In_state = In_state
        self.Out_state = Out_state
        self.Pressure_drop = Pressure_drop
        self.Temperature_drop = Temperature_drop
        super().__init__(Model, ID, Calculate)

    @validate_inputs
    def Cal(self):
        self.In = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]
        
        # Resolve mass flow rate
        mass, deffult = self._resolve_mass_flowrate(self.In, self.Out)
        self.In.Mass_flowrate = self.Model.Point[self.In_state].Mass_flowrate = mass
        self.Out.Mass_flowrate = self.Model.Point[self.Out_state].Mass_flowrate = mass

        # ===== PRESSURE SOLVER =====
        # Need at least 2 of: P_in, P_out, ΔP
        known_P_count = sum([
            self.In.P is not None,
            self.Out.P is not None,
            self.Pressure_drop is not None
        ])
        
        if known_P_count < 2:
            raise ValueError(
                f"Pipe {self.ID}: Need at least 2 of (P_in, P_out, Pressure_drop)"
            )
        
        if self.In.P is None:
            self.In.P = self.Out.P + self.Pressure_drop
        elif self.Out.P is None:
            self.Out.P = self.In.P - self.Pressure_drop
        else:  # Both pressures known
            calculated_drop = self.In.P - self.Out.P
            if self.Pressure_drop is None:
                self.Pressure_drop = calculated_drop
            elif abs(calculated_drop - self.Pressure_drop) > self.TOLERANCE_PRESSURE:
                raise ValueError(
                    f"Pipe {self.ID}: Pressure inconsistency. "
                    f"Specified drop: {self.Pressure_drop/1e5:.3f} bar, "
                    f"Calculated: {calculated_drop/1e5:.3f} bar"
                )
        
        # Validate positive pressures
        if self.In.P < 0 or self.Out.P < 0:
            raise ValueError(
                f"Pipe {self.ID}: Negative pressure encountered "
                f"(P_in={self.In.P/1e5:.2f}, P_out={self.Out.P/1e5:.2f} bar)"
            )

        # ===== TEMPERATURE SOLVER =====
        # Need at least 2 of: T_in, T_out, ΔT
        known_T_count = sum([
            self.In.T is not None,
            self.Out.T is not None,
            self.Temperature_drop is not None
        ])
        
        if known_T_count < 2:
            raise ValueError(
                f"Pipe {self.ID}: Need at least 2 of (T_in, T_out, Temperature_drop)"
            )
        
        if self.In.T is None:
            self.In.T = self.Out.T + self.Temperature_drop
        elif self.Out.T is None:
            self.Out.T = self.In.T - self.Temperature_drop
        else:  # Both temperatures known
            calculated_drop = self.In.T - self.Out.T
            if self.Temperature_drop is None:
                self.Temperature_drop = calculated_drop
            elif abs(calculated_drop - self.Temperature_drop) > self.TOLERANCE_TEMPERATURE:
                raise ValueError(
                    f"Pipe {self.ID}: Temperature inconsistency. "
                    f"Specified drop: {self.Temperature_drop:.2f} K, "
                    f"Calculated: {calculated_drop:.2f} K"
                )

        # Calculate exergy destruction
        try:
            self.Ex_D = self.In.Ex - self.Out.Ex
            if self.Ex_D < 0:
                warnings.warn(
                    f"Pipe {self.ID}: Negative exergy destruction ({self.Ex_D:.2f} W). "
                    f"Check dead state definition."
                )
        except (AttributeError, TypeError):
            self.Ex_D = None

        self.Solution_Status = True
        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Pipe): not yet solved."
        return (
            f"{self.ID} (Pipe):\n"
            f"  P_in  → P_out  : {self.In.P/1e5:.2f} → {self.Out.P/1e5:.2f} bar "
            f"(ΔP = {self.Pressure_drop/1e5:.3f} bar)\n"
            f"  T_in  → T_out  : {self.In.T:.2f} → {self.Out.T:.2f} K "
            f"(ΔT = {self.Temperature_drop:.2f} K)\n"
            f"  Mass flow      : {self.In.Mass_flowrate:.4f} kg/s\n"
            f"  Exergy destr.  : {self.Ex_D}\n"
            f"  Solved         : {self.Solution_Status}"
        )


# ================================================================== #
#   EXPANSION VALVE
# ================================================================== #
class Expansion_valve(Component):
    """Isenthalpic expansion device (throttling valve)"""
    
    TOLERANCE_ENTHALPY = 1.0  # J/kg

    def __init__(self, Model, ID, In_state, Out_state, Calculate=False):
        self.In_state = In_state
        self.Out_state = Out_state
        super().__init__(Model, ID, Calculate)

    @validate_inputs
    def Cal(self):
        self.In = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]
        
        # Resolve mass flow rate
        mass, default = self._resolve_mass_flowrate(self.In, self.Out)
        self.In.Mass_flowrate = self.Model.Point[self.In_state].Mass_flowrate = mass
        self.Out.Mass_flowrate = self.Model.Point[self.Out_state].Mass_flowrate = mass
        # Validate pressures exist
        if self.In.P is None and self.Out.P is None:
            raise ValueError(
                f"Expansion_valve {self.ID}: At least one pressure must be specified"
            )
        
        if self.In.P is not None and self.Out.P is not None:
            if self.In.P <= self.Out.P:
                warnings.warn(
                    f"Expansion_valve {self.ID}: Outlet pressure ({self.Out.P/1e5:.2f} bar) "
                    f"is not less than inlet pressure ({self.In.P/1e5:.2f} bar)"
                )

        # ===== ISENTHALPIC PROCESS (H_in = H_out) =====
        if self.In.H is None and self.Out.H is None:
            raise ValueError(
                f"Expansion_valve {self.ID}: At least one enthalpy must be known"
            )
        
        if self.In.H is None:
            # Calculate inlet from outlet
            if self.In.P is None:
                raise ValueError(
                    f"Expansion_valve {self.ID}: Cannot determine inlet state. "
                    f"Need P_in when H_in is unknown"
                )
            self.In = self.Model.Prop(
                self.In.fluid,
                StatePointName=self.In.StatePointName,
                P=self.In.P, H=self.Out.H,
                Mass_flowrate=mass
            )
            
        elif self.Out.H is None:
            # Calculate outlet from inlet
            if self.Out.P is None:
                raise ValueError(
                    f"Expansion_valve {self.ID}: Cannot determine outlet state. "
                    f"Need P_out when H_out is unknown"
                )
            self.Out = self.Model.Prop(
                self.Out.fluid,
                StatePointName=self.Out.StatePointName,
                P=self.Out.P, H=self.In.H,
                Mass_flowrate=mass
            )
            
        else:
            # Both enthalpies specified - check consistency
            if abs(self.In.H - self.Out.H) > self.TOLERANCE_ENTHALPY:
                raise ValueError(
                    f"Expansion_valve {self.ID}: Enthalpy mismatch in isenthalpic process. "
                    f"H_in = {self.In.H:.2f} J/kg, H_out = {self.Out.H:.2f} J/kg, "
                    f"ΔH = {abs(self.In.H - self.Out.H):.2f} J/kg"
                )

        # Calculate exergy destruction (should be positive for irreversible throttling)
        try:
            self.Ex_D = self.In.Ex - self.Out.Ex
            if self.Ex_D < 0:
                warnings.warn(
                    f"Expansion_valve {self.ID}: Negative exergy destruction ({self.Ex_D:.2f} W)"
                )
        except (AttributeError, TypeError):
            self.Ex_D = None

        self.Solution_Status = True
        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Expansion Valve): not yet solved."
        h_avg = ((self.In.H + self.Out.H) / 2
                 if None not in (self.In.H, self.Out.H) else None)
        return (
            f"{self.ID} (Expansion Valve):\n"
            f"  P_in  → P_out  : {self.In.P/1e5:.2f} → {self.Out.P/1e5:.2f} bar "
            f"(ΔP = {(self.In.P - self.Out.P)/1e5:.3f} bar)\n"
            f"  Enthalpy       : {h_avg:.2f} J/kg (isenthalpic)\n"
            f"  Mass flow      : {self.In.Mass_flowrate:.4f} kg/s\n"
            f"  Exergy destr.  : {self.Ex_D}\n"
            f"  Solved         : {self.Solution_Status}"
        )


# ================================================================== #
#   COMPONENT TOLERANCE CONFIGURATION
# ================================================================== #
class ComponentConfig:
    """Global configuration for component tolerances"""
    
    @staticmethod
    def set_tolerances(**kwargs):
        """
        Set tolerances for all components.
        
        Example:
            ComponentConfig.set_tolerances(
                pressure=1e4,           # 10 kPa
                temperature=0.5,        # 0.5 K
                enthalpy=10.0           # 10 J/kg
            )
        """
        mapping = {
            'pressure': ['Pipe'],
            'temperature': ['Pipe'],
            'enthalpy': ['Mixer', 'Expansion_valve'],
            'fraction_sum': ['Splitter']
        }
        
        for param, value in kwargs.items():
            if param in mapping:
                for component_name in mapping[param]:
                    component_class = globals()[component_name]
                    attr_name = f'TOLERANCE_{param.upper()}'
                    if hasattr(component_class, attr_name):
                        setattr(component_class, attr_name, value)
                        print(f"Set {component_name}.{attr_name} = {value}")
            else:
                warnings.warn(f"Unknown tolerance parameter: {param}")