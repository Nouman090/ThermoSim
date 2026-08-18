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
                 n_isen=1.0, n_mech=1.0, Calculate=False):
        self.In_state  = In_state
        self.Out_state = Out_state
        self.n_isen = n_isen
        self.n_mech = n_mech
        self.work   = 0.0
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]

        # Resolve the mass flow for this component.  A straight-through
        # machine has one stream, so inlet and outlet always carry the same
        # flow — we do NOT need the model-wide propagation graph to know that.
        # Using `or` here would treat a legitimate 0.0 kg/s as missing, so the
        # test is explicit.
        mass, used_default = self._resolve_mass_flowrate(self.In, self.Out)

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

        _propagate_mass(self, mass, used_default)

        w_specific = self.In.H - self.Out.H
        self.work  = (w_specific * mass * self.n_mech) if mass is not None else None

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
                 n_isen=1.0, n_mech=1.0, Calculate=False):
        self.In_state  = In_state
        self.Out_state = Out_state
        self.n_isen = n_isen
        self.n_mech = n_mech
        self.work   = 0.0
        self.pressure_ratio = None
        super().__init__(Model, ID, Calculate)

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]

        # Resolve the mass flow for this component.  A straight-through
        # machine has one stream, so inlet and outlet always carry the same
        # flow — we do NOT need the model-wide propagation graph to know that.
        # Using `or` here would treat a legitimate 0.0 kg/s as missing, so the
        # test is explicit.
        mass, used_default = self._resolve_mass_flowrate(self.In, self.Out)

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

        _propagate_mass(self, mass, used_default)

        w_specific = self.Out.H - self.In.H
        self.work = (w_specific * mass / self.n_mech) if mass is not None else None
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
                 n_isen=1.0, n_mech=1.0, Calculate=False):
        self.In_state  = In_state
        self.Out_state = Out_state
        self.n_isen = n_isen
        self.n_mech = n_mech
        self.work   = 0.0
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

        # Resolve the mass flow for this component.  A straight-through
        # machine has one stream, so inlet and outlet always carry the same
        # flow — we do NOT need the model-wide propagation graph to know that.
        # Using `or` here would treat a legitimate 0.0 kg/s as missing, so the
        # test is explicit.
        mass, used_default = self._resolve_mass_flowrate(self.In, self.Out)

        if self.Compressibility == 'Compressible':
            w_specific = self._solve_compressible()
        elif self.Compressibility == 'Incompressible':
            w_specific = self._solve_incompressible()
        else:
            raise ValueError(
                f"Invalid Compressibility '{self.Compressibility}' in {self.ID}. "
                "Use 'Compressible' or 'Incompressible'."
            )

        _propagate_mass(self, mass, used_default)

        try:
            self.work = (w_specific * mass / self.n_mech) if mass is not None else None
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
