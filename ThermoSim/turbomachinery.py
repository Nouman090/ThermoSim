"""
turbomachinery.py
-----------------
Turbine, Compressor, and Pump components.


"""

import warnings
from scipy.optimize import fsolve
from .base_component import Component


# ------------------------------------------------------------------ #
#  Shared helpers
# ------------------------------------------------------------------ #

def _check_convergence(sol, info, msg_prefix):
    """Raise a warning if fsolve did not converge (ier != 1)."""
    ier = info[4] if isinstance(info, tuple) else 1
    if ier != 1:
        warnings.warn(
            f"{msg_prefix}: fsolve did not converge (ier={ier}). "
            "Result may be unreliable. Try a better initial guess.",
            RuntimeWarning,
            stacklevel=3,
        )


def _resolve_mass_with_target_work(component, In, Out, work_given, w_eff):
    """
    Resolve a turbomachine's mass flow rate, accounting for an optional
    target ``work`` (the constructor's ``work=`` argument) alongside the
    usual port-based resolution.

    ``w_eff`` is the mechanical-efficiency-adjusted specific work implied by
    the already-solved (In.H, Out.H) pair -- turbine: (h_in-h_out)*n_mech;
    pump/compressor: (h_out-h_in)/n_mech -- i.e. exactly what multiplying by
    mass gives ``self.work``. It must be computed AFTER the enthalpy solve,
    which is why this runs later in Cal() than the old direct call to
    ``_resolve_mass_flowrate`` used to.

    Two cases, matching ``base_component._resolve_mass_flowrate``'s
    ``(mass, used_default)`` contract:

      * Neither port carries a mass flow rate, and ``work_given`` is not
        None -> mass is DERIVED from ``work_given / w_eff`` instead of
        falling back to the 1 kg/s guess (and without that fallback's
        warning -- a real value was available, just not on a port).
      * A mass flow rate IS already resolvable from the ports -> that value
        is authoritative; ``work_given``, if supplied, is only checked for
        consistency against the work it implies, exactly like the existing
        "both In.H and Out.H known" branches check enthalpy consistency.
    """
    both_missing = In.Mass_flowrate is None and Out.Mass_flowrate is None

    if work_given is not None and both_missing:
        if abs(w_eff) < 1e-9:
            raise ValueError(
                f"{component.ID}: cannot derive a mass flow rate from "
                f"work={work_given:.2f} W -- the solved state pair implies "
                f"~zero specific work."
            )
        mass = work_given / w_eff
        if mass <= 0:
            raise ValueError(
                f"{component.ID}: work={work_given:.2f} W implies a "
                f"non-positive mass flow rate ({mass:.6g} kg/s) for the "
                f"solved state pair. Check the sign of `work` against the "
                f"device (turbine work is positive OUT; pump/compressor "
                f"work is positive IN)."
            )
        return mass, False

    mass, used_default = component._resolve_mass_flowrate(In, Out)

    if work_given is not None and not used_default:
        expected = mass * w_eff
        rel_err = abs(work_given - expected) / max(abs(expected), 1.0)
        if rel_err > 1e-3:
            warnings.warn(
                f"{component.ID}: supplied work={work_given:.2f} W is "
                f"inconsistent with the mass flow resolved from the state "
                f"points ({mass:.6g} kg/s implies {expected:.2f} W, "
                f"{rel_err:.2%} difference). The port-resolved mass flow "
                f"is used as-is; `work` is not re-applied.",
                UserWarning, stacklevel=3
            )

    return mass, used_default


def _propagate_mass(component, mass, used_default=False):
    """
    Write the resolved mass flow onto both ports of a single-stream machine.

    Turbine, Compressor and Pump each carry exactly one stream, so the inlet
    and outlet mass flows are equal by definition.  Up to v3.2.1 the flow was
    kept as a local variable and never written back, on the assumption that
    ``ThermodynamicModel``'s flow-propagation graph would fill it in.  But
    that graph is only built inside ``Solve()`` / ``enable_flow_propagation()``,
    so the documented ``Calculate=True`` workflow left every downstream point
    at ``Mass_flowrate=None`` and the next heat exchanger could not solve.

    Writing it here is safe alongside the graph: assignment goes through
    ``Prop.__setattr__``, which fires the propagator when one is attached, so
    branch topology is still handled by the network.

    ``used_default`` must be forwarded from ``_resolve_mass_flowrate``.  A
    defaulted flow is a GUESS (the 1 kg/s fallback used when no port carries a
    rate) and must never be written to a shared state point: a later component
    that resolves the real rate would then see its own correct value disagree
    with the guess and raise a spurious mismatch.  Work for this component is
    still evaluated on the fallback, which is the documented lenient
    behaviour -- only the propagation is suppressed.
    """
    if mass is None or used_default:
        return
    for port in (component.In, component.Out):
        if port is not None and port.Mass_flowrate is None:
            port.Mass_flowrate = mass


def _isentropic_work_turbine(h_in, h_isen, n_isen):
    """Actual specific work output of a turbine (positive value)."""
    return (h_in - h_isen) * n_isen


def _isentropic_work_pump_compressor(h_in, h_isen, n_isen):
    """Actual specific work input of a pump/compressor (positive value)."""
    return (h_isen - h_in) / n_isen


def _warn_consistency(component_id, h_in, h_out, h_isen, n_isen, kind):
    """
    Warn if the supplied (h_in, h_out) pair is inconsistent with n_isen
    for a turbine or a pump/compressor.

    kind : 'turbine' | 'pump'
    """
    if kind == 'turbine':
        expected_h_out = h_in - (h_in - h_isen) * n_isen
    else:  # pump / compressor
        expected_h_out = h_in + (h_isen - h_in) / n_isen

    rel_err = abs(h_out - expected_h_out) / (abs(expected_h_out) + 1e-12)
    if rel_err > 1e-4:
        warnings.warn(
            f"{component_id}: supplied h_in={h_in:.2f} and h_out={h_out:.2f} "
            f"are inconsistent with η_isen={n_isen:.4f} "
            f"(expected h_out≈{expected_h_out:.2f}, rel. error={rel_err:.2%}). "
            "The isentropic efficiency is NOT re-applied; existing enthalpies "
            "are used as-is.",
            UserWarning,
            stacklevel=3,
        )


# ================================================================== #
#   TURBINE
# ================================================================== #
class Turbine(Component):

    def __init__(self, Model, ID, In_state, Out_state,
                 n_isen=1.0, n_mech=1.0, work=None, Calculate=False):
        self.In_state  = In_state
        self.Out_state = Out_state
        self.n_isen = n_isen
        self.n_mech = n_mech
        self.work   = 0.0
        # The constructor's target work, if any -- kept separate from
        # self.work (which Cal() always overwrites with the solved output)
        # so a Solve() round-robin can check it on every call regardless of
        # what the previous iteration computed.
        self._work_given = work
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]

        if self.Out.H is None and self.In.H is not None:
            # ---- forward: known inlet, solve outlet ----
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_isen_tmp'
            )
            h_out = self.In.H - _isentropic_work_turbine(
                self.In.H, h_isen_pt.H, self.n_isen)
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=h_out,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )

        elif self.In.H is None and self.Out.H is not None:
            # ---- reverse: known outlet, solve inlet ----
            h_guess = self.Out.H / self.n_isen

            def _eq(h_in):
                _in = self.Model.Prop(
                    self.In.fluid, P=self.In.P, H=h_in[0],
                    StatePointName=self.In.StatePointName
                )
                h_is = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, S=_in.S,
                    StatePointName='_isen_tmp'
                )
                return [self.Out.H - h_in[0] +
                        _isentropic_work_turbine(h_in[0], h_is.H, self.n_isen)]

            sol, info, ier, msg = fsolve(_eq, [h_guess], full_output=True)
            _check_convergence(sol, (None, None, None, None, ier),
                               f"{self.ID} reverse solve")
            h_in_sol = sol[0]

            if h_in_sol <= self.Out.H:
                warnings.warn(
                    f"{self.ID}: reverse solve yielded h_in ({h_in_sol:.2f}) ≤ "
                    f"h_out ({self.Out.H:.2f}). Check inlet pressure and η_isen.",
                    UserWarning, stacklevel=2
                )

            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=h_in_sol,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )

        elif self.In.H is not None and self.Out.H is not None:
            # ---- both known: rebuild full property sets + consistency check ----
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_isen_tmp'
            )
            _warn_consistency(self.ID, self.In.H, self.Out.H,
                              h_isen_pt.H, self.n_isen, 'turbine')
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=self.In.H,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=self.Out.H,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )

        else:
            raise ValueError(
                f"{self.ID}: cannot solve — both In.H and Out.H are None. "
                "Provide at least one enthalpy (or a full state) to proceed."
            )

        w_specific = self.In.H - self.Out.H
        w_eff = w_specific * self.n_mech
        mass, used_default = _resolve_mass_with_target_work(
            self, self.In, self.Out, self._work_given, w_eff)

        _propagate_mass(self, mass, used_default)

        self.work  = (w_eff * mass) if mass is not None else None

        try:
            self.Ex_D = self.In.Ex - self.Out.Ex - self.work
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True

        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Turbine): not yet solved."
        try:
            return (
                f"{self.ID} (Turbine):\n"
                f"  P_in            : {self.In.P/1e5:.2f} bar\n"
                f"  P_out           : {self.Out.P/1e5:.2f} bar\n"
                f"  η_isen          : {self.n_isen*100:.1f} %\n"
                f"  η_mech          : {self.n_mech*100:.1f} %\n"
                f"  Work            : {self.work:.2f} W\n"
                f"  Exergy dest.    : {self.Ex_D}\n"
                f"  Solved          : {self.Solution_Status}"
            )
        except:
            return (
                f"{self.ID} (Turbine):\n"
                f"  P_in            : {self.In.P} bar\n"
                f"  P_out           : {self.Out.P} bar\n"
                f"  η_isen          : {self.n_isen*100:.1f} %\n"
                f"  η_mech          : {self.n_mech*100:.1f} %\n"
                f"  Work            : {self.work} W\n"
                f"  Exergy dest.    : {self.Ex_D}\n"
                f"  Solved          : {self.Solution_Status}"
            )


# ================================================================== #
#   COMPRESSOR  (for gas cycles, e.g. Brayton)
# ================================================================== #
class Compressor(Component):
    """
    Adiabatic gas compressor with isentropic efficiency.

    Identical physics to Pump but intended for gas-phase working fluids
    where density changes significantly.  Also tracks pressure ratio.

    Parameters
    ----------
    Model : ThermodynamicModel
    ID : str
    In_state, Out_state : str
        State-point names.
    n_isen : float
        Isentropic efficiency (0-1).
    n_mech : float
        Mechanical efficiency (0-1).
    Calculate : bool
        If True, solve immediately upon creation.
    """

    def __init__(self, Model, ID, In_state, Out_state,
                 n_isen=1.0, n_mech=1.0, work=None, Calculate=False):
        self.In_state  = In_state
        self.Out_state = Out_state
        self.n_isen = n_isen
        self.n_mech = n_mech
        self.work   = 0.0
        self._work_given = work
        self.pressure_ratio = None
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]

        if self.Out.H is None and self.In.H is not None:
            # ---- forward: known inlet, solve outlet ----
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_comp_isen_tmp'
            )
            h_out = self.In.H + _isentropic_work_pump_compressor(
                self.In.H, h_isen_pt.H, self.n_isen)
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=h_out,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )

        elif self.In.H is None and self.Out.H is not None:
            # ---- reverse: known outlet, solve inlet ----
            h_guess = self.Out.H * self.n_isen

            def _eq(h_in):
                _in = self.Model.Prop(
                    self.In.fluid, P=self.In.P, H=h_in[0],
                    StatePointName=self.In.StatePointName
                )
                h_is = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, S=_in.S,
                    StatePointName='_comp_isen_tmp'
                )
                return [self.Out.H - h_in[0] -
                        _isentropic_work_pump_compressor(h_in[0], h_is.H, self.n_isen)]

            sol, info, ier, msg = fsolve(_eq, [h_guess], full_output=True)
            _check_convergence(sol, (None, None, None, None, ier),
                               f"{self.ID} reverse solve")
            h_in_sol = sol[0]

            if h_in_sol >= self.Out.H:
                warnings.warn(
                    f"{self.ID}: reverse solve yielded h_in ({h_in_sol:.2f}) ≥ "
                    f"h_out ({self.Out.H:.2f}). Check outlet pressure and η_isen.",
                    UserWarning, stacklevel=2
                )

            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=h_in_sol,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )

        elif self.In.H is not None and self.Out.H is not None:
            # ---- both known: rebuild + consistency check ----
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_comp_isen_tmp'
            )
            _warn_consistency(self.ID, self.In.H, self.Out.H,
                              h_isen_pt.H, self.n_isen, 'pump')
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=self.In.H,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=self.Out.H,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )

        else:
            raise ValueError(
                f"{self.ID}: cannot solve — both In.H and Out.H are None."
            )

        w_specific = self.Out.H - self.In.H
        w_eff = w_specific / self.n_mech
        mass, used_default = _resolve_mass_with_target_work(
            self, self.In, self.Out, self._work_given, w_eff)

        _propagate_mass(self, mass, used_default)

        self.work = (w_eff * mass) if mass is not None else None
        self.pressure_ratio = self.Out.P / self.In.P

        try:
            self.Ex_D = self.work + self.In.Ex - self.Out.Ex
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True
        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Compressor): not yet solved."
        try:
            return (
                f"{self.ID} (Compressor):\n"
                f"  P_in            : {self.In.P/1e5:.2f} bar\n"
                f"  P_out           : {self.Out.P/1e5:.2f} bar\n"
                f"  η_isen          : {self.n_isen*100:.1f} %\n"
                f"  η_mech          : {self.n_mech*100:.1f} %\n"
                f"  Work            : {self.work:.2f} W\n"
                f"  Exergy dest.    : {self.Ex_D}\n"
                f"  Solved          : {self.Solution_Status}"
            )
        except:
            return (
                f"{self.ID} (Compressor):\n"
                f"  P_in            : {self.In.P} bar\n"
                f"  P_out           : {self.Out.P} bar\n"
                f"  η_isen          : {self.n_isen*100:.1f} %\n"
                f"  η_mech          : {self.n_mech*100:.1f} %\n"
                f"  Work            : {self.work} W\n"
                f"  Exergy dest.    : {self.Ex_D}\n"
                f"  Solved          : {self.Solution_Status}"
            )


# ================================================================== #
#   PUMP
# ================================================================== #
class Pump(Component):

    def __init__(self, Model, ID, In_state, Out_state,
                 Compressibility='Compressible',
                 n_isen=1.0, n_mech=1.0, work=None, Calculate=False):
        self.In_state  = In_state
        self.Out_state = Out_state
        self.n_isen = n_isen
        self.n_mech = n_mech
        self.work   = 0.0
        self._work_given = work
        self.Compressibility = Compressibility
        super().__init__(Model, ID, Calculate)

    # ---------------------------------------------------------------- #
    #  Internal sub-solvers — keeps Cal() readable
    # ---------------------------------------------------------------- #

    def _solve_compressible(self):
        """Solve compressible pump in all three enthalpy configurations."""
        if self.Out.H is None and self.In.H is not None:
            # forward
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_isen_tmp', Mass_flowrate=self.In.Mass_flowrate
            )
            h_out = self.In.H + _isentropic_work_pump_compressor(
                self.In.H, h_isen_pt.H, self.n_isen)
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=h_out,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )
            return self.Out.H - self.In.H  # w_specific

        elif self.In.H is None and self.Out.H is not None:
            # reverse
            h_guess = self.Out.H * self.n_isen

            def _eq(h_in):
                _in = self.Model.Prop(
                    self.In.fluid, P=self.In.P, H=h_in[0],
                    StatePointName=self.In.StatePointName
                )
                h_is = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, S=_in.S,
                    StatePointName='_isen_tmp'
                )
                return [self.Out.H - h_in[0] -
                        _isentropic_work_pump_compressor(h_in[0], h_is.H, self.n_isen)]

            sol, info, ier, msg = fsolve(_eq, [h_guess], full_output=True)
            _check_convergence(sol, (None, None, None, None, ier),
                               f"{self.ID} reverse solve")
            h_in_sol = sol[0]

            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=h_in_sol,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )
            return self.Out.H - self.In.H

        elif self.In.H is not None and self.Out.H is not None:
            # both known — consistency check then rebuild
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_isen_tmp'
            )
            _warn_consistency(self.ID, self.In.H, self.Out.H,
                              h_isen_pt.H, self.n_isen, 'pump')
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=self.In.H,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=self.Out.H,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )
            return self.Out.H - self.In.H

        else:
            raise ValueError(
                f"{self.ID} (Compressible): both In.H and Out.H are None."
            )

    def _solve_incompressible(self):
        """
        Incompressible pump — density assumed constant at inlet value.

        Handles three cases:
          forward  — In.H known, Out.H unknown
          reverse  — Out.H known, In.H unknown  (NEW)
          both     — both known, consistency check  (NEW)
        """
        if self.In.D is None:
            raise ValueError(
                f"{self.ID} (Incompressible): inlet density (In.D) is None. "
                "The inlet state must be fully defined before calling Cal()."
            )
        D = self.In.D

        if self.Out.H is None and self.In.H is not None:
            # forward
            w_ideal   = (self.Out.P - self.In.P) / D
            w_specific = w_ideal / self.n_isen
            h_out = self.In.H + w_specific
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=h_out,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )
            return w_specific

        elif self.In.H is None and self.Out.H is not None:
            # reverse — compute inlet enthalpy from outlet
            w_ideal   = (self.Out.P - self.In.P) / D
            w_specific = w_ideal / self.n_isen
            h_in = self.Out.H - w_specific
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=h_in,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )
            return w_specific

        elif self.In.H is not None and self.Out.H is not None:
            # both known — check consistency
            w_ideal    = (self.Out.P - self.In.P) / D
            w_specific = w_ideal / self.n_isen
            expected_h_out = self.In.H + w_specific
            rel_err = abs(self.Out.H - expected_h_out) / (abs(expected_h_out) + 1e-12)
            if rel_err > 1e-4:
                warnings.warn(
                    f"{self.ID} (Incompressible): supplied h_out={self.Out.H:.2f} "
                    f"inconsistent with incompressible model (expected "
                    f"h_out≈{expected_h_out:.2f}, rel. error={rel_err:.2%}).",
                    UserWarning, stacklevel=3
                )
            # Rebuild both points
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=self.In.H,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=self.In.Mass_flowrate
            )
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=self.Out.H,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=self.Out.Mass_flowrate
            )
            return self.Out.H - self.In.H

        else:
            raise ValueError(
                f"{self.ID} (Incompressible): both In.H and Out.H are None."
            )

    # ---------------------------------------------------------------- #
    #  Main calculation entry point
    # ---------------------------------------------------------------- #

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]

        if self.Compressibility == 'Compressible':
            w_specific = self._solve_compressible()
        elif self.Compressibility == 'Incompressible':
            w_specific = self._solve_incompressible()
        else:
            raise ValueError(
                f"Invalid Compressibility '{self.Compressibility}' in {self.ID}. "
                "Use 'Compressible' or 'Incompressible'."
            )

        w_eff = w_specific / self.n_mech
        mass, used_default = _resolve_mass_with_target_work(
            self, self.In, self.Out, self._work_given, w_eff)

        _propagate_mass(self, mass, used_default)

        try:
            self.work = (w_eff * mass) if mass is not None else None
            self.Ex_D = self.work + self.In.Ex - self.Out.Ex
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True
        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        if not getattr(self, 'Solution_Status', False):
            return f"{self.ID} (Pump): not yet solved."
        try:
            return (
                f"{self.ID} (Pump):\n"
                f"  Compressibility : {self.Compressibility}\n"
                f"  P_in            : {self.In.P/1e5:.2f} bar\n"
                f"  P_out           : {self.Out.P/1e5:.2f} bar\n"
                f"  η_isen          : {self.n_isen*100:.1f} %\n"
                f"  η_mech          : {self.n_mech*100:.1f} %\n"
                f"  Work            : {self.work:.2f} W\n"
                f"  Exergy dest.    : {self.Ex_D}\n"
                f"  Solved          : {self.Solution_Status}"
            )
        except:
            return (
                f"{self.ID} (Pump):\n"
                f"  Compressibility : {self.Compressibility}\n"
                f"  P_in            : {self.In.P} bar\n"
                f"  P_out           : {self.Out.P} bar\n"
                f"  η_isen          : {self.n_isen*100:.1f} %\n"
                f"  η_mech          : {self.n_mech*100:.1f} %\n"
                f"  Work            : {self.work} W\n"
                f"  Exergy dest.    : {self.Ex_D}\n"
                f"  Solved          : {self.Solution_Status}"
            )
