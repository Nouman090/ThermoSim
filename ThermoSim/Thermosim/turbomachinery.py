"""
turbomachinery.py
-----------------
Turbine and Pump components.
"""

import warnings
from scipy.optimize import fsolve
from .base_component import Component


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

        mass, used_default = self._resolve_mass_flowrate(self.In, self.Out)

        if self.Out.H is None:
            # --- forward: known inlet, unknown outlet ---
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_isen_tmp'
            )
            h_out = self.In.H - (self.In.H - h_isen_pt.H) * self.n_isen
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=h_out,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=mass
            )

        elif self.In.H is None:
            # --- reverse: known outlet, unknown inlet ---
            def _eq(h_in):
                _in = self.Model.Prop(
                    self.In.fluid, P=self.In.P, H=h_in[0],
                    StatePointName=self.In.StatePointName
                )
                h_is = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, S=_in.S,
                    StatePointName='_isen_tmp'
                )
                return self.Out.H - h_in[0] + (h_in[0] - h_is.H) * self.n_isen

            h_in = fsolve(_eq, [self.Out.H])
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=h_in[0],
                StatePointName=self.In.StatePointName,
                Mass_flowrate=mass
            )

        else:
            # --- both known: just rebuild full property set ---
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=self.In.H,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=mass
            )
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=self.Out.H,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=mass
            )

        w_specific = self.In.H - self.Out.H
        self.work  = w_specific * mass * self.n_mech

        # exergy destruction
        try:
            self.Ex_D = self.In.Ex - self.Out.Ex - self.work
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True

        if used_default:
            self.In.Mass_flowrate  = None
            self.Out.Mass_flowrate = None

        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        return (
            f"{self.ID} (Turbine):\n"
            f"  P_in        : {self.In.P/1e5:.2f} bar\n"
            f"  P_out       : {self.Out.P/1e5:.2f} bar\n"
            f"  η_isen      : {self.n_isen*100:.1f} %\n"
            f"  η_mech      : {self.n_mech*100:.1f} %\n"
            f"  Work        : {self.work:.2f} W\n"
            f"  Exergy dest.: {self.Ex_D}\n"
            f"  Solved      : {self.Solution_Status}"
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

        mass, used_default = self._resolve_mass_flowrate(self.In, self.Out)

        if self.Out.H is None:
            # --- forward: known inlet, unknown outlet ---
            h_isen_pt = self.Model.Prop(
                self.In.fluid, P=self.Out.P, S=self.In.S,
                StatePointName='_comp_isen_tmp'
            )
            h_out = self.In.H + (h_isen_pt.H - self.In.H) / self.n_isen
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=h_out,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=mass
            )

        elif self.In.H is None:
            # --- reverse: known outlet, unknown inlet ---
            from scipy.optimize import fsolve

            def _eq(h_in):
                _in = self.Model.Prop(
                    self.In.fluid, P=self.In.P, H=h_in[0],
                    StatePointName=self.In.StatePointName
                )
                h_is = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, S=_in.S,
                    StatePointName='_comp_isen_tmp'
                )
                return self.Out.H - h_in[0] - (h_is.H - h_in[0]) / self.n_isen

            h_in = fsolve(_eq, [self.Out.H])
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=h_in[0],
                StatePointName=self.In.StatePointName,
                Mass_flowrate=mass
            )

        else:
            # --- both known ---
            self.In = self.Model.Prop(
                self.In.fluid, P=self.In.P, H=self.In.H,
                StatePointName=self.In.StatePointName,
                Mass_flowrate=mass
            )
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=self.Out.H,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=mass
            )

        w_specific = self.Out.H - self.In.H
        self.work = w_specific * mass / self.n_mech
        self.pressure_ratio = self.Out.P / self.In.P

        # exergy: work input, so Ex_D = W_in + Ex_in - Ex_out
        try:
            self.Ex_D = self.work + self.In.Ex - self.Out.Ex
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True

        if used_default:
            self.In.Mass_flowrate  = None
            self.Out.Mass_flowrate = None

        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        return (
            f"{self.ID} (Compressor):\n"
            f"  P_in             : {self.In.P/1e5:.2f} bar\n"
            f"  P_out            : {self.Out.P/1e5:.2f} bar\n"
            f"  Pressure ratio   : {self.pressure_ratio:.2f}\n"
            f"  η_isen           : {self.n_isen*100:.1f} %\n"
            f"  η_mech           : {self.n_mech*100:.1f} %\n"
            f"  Work input       : {self.work:.2f} W  ({self.work/1e3:.2f} kW)\n"
            f"  Exergy destr.    : {self.Ex_D}\n"
            f"  Solved           : {self.Solution_Status}"
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

    def Cal(self):
        self.In  = self.Model.Point[self.In_state]
        self.Out = self.Model.Point[self.Out_state]

        mass, used_default = self._resolve_mass_flowrate(self.In, self.Out)

        if self.Compressibility == 'Compressible':
            if self.Out.H is None:
                h_isen_pt = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, S=self.In.S,
                    StatePointName='_isen_tmp', Mass_flowrate=mass
                )
                h_out = self.In.H + (h_isen_pt.H - self.In.H) / self.n_isen
                self.Out = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, H=h_out,
                    StatePointName=self.Out.StatePointName,
                    Mass_flowrate=mass
                )

            elif self.In.H is None:
                def _eq(h_in):
                    _in = self.Model.Prop(
                        self.In.fluid, P=self.In.P, H=h_in[0],
                        StatePointName=self.In.StatePointName
                    )
                    h_is = self.Model.Prop(
                        self.In.fluid, P=self.Out.P, S=_in.S,
                        StatePointName='_isen_tmp'
                    )
                    return self.Out.H - h_in[0] - (h_is.H - h_in[0]) / self.n_isen

                h_in = fsolve(_eq, [self.Out.H])
                self.In = self.Model.Prop(
                    self.In.fluid, P=self.In.P, H=h_in[0],
                    StatePointName=self.In.StatePointName,
                    Mass_flowrate=mass
                )
            else:
                self.In = self.Model.Prop(
                    self.In.fluid, P=self.In.P, H=self.In.H,
                    StatePointName=self.In.StatePointName,
                    Mass_flowrate=mass
                )
                self.Out = self.Model.Prop(
                    self.In.fluid, P=self.Out.P, H=self.Out.H,
                    StatePointName=self.Out.StatePointName,
                    Mass_flowrate=mass
                )

            w_specific = self.Out.H - self.In.H

        elif self.Compressibility == 'Incompressible':
            D = self.In.D
            w_ideal = (self.Out.P - self.In.P) / D
            w_specific = w_ideal / self.n_isen
            self.Out.H = self.In.H + w_specific
            self.Out = self.Model.Prop(
                self.In.fluid, P=self.Out.P, H=self.Out.H,
                StatePointName=self.Out.StatePointName,
                Mass_flowrate=mass
            )
        else:
            raise ValueError(
                f"Invalid Compressibility '{self.Compressibility}' in {self.ID}. "
                f"Use 'Compressible' or 'Incompressible'."
            )

        self.work = w_specific * mass / self.n_mech

        try:
            self.Ex_D = self.work + self.In.Ex - self.Out.Ex
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True

        if used_default:
            self.In.Mass_flowrate  = None
            self.Out.Mass_flowrate = None

        self._update_model_points(In_state=self.In, Out_state=self.Out)

    def __str__(self):
        return (
            f"{self.ID} (Pump):\n"
            f"  P_in        : {self.In.P/1e5:.2f} bar\n"
            f"  P_out       : {self.Out.P/1e5:.2f} bar\n"
            f"  η_isen      : {self.n_isen*100:.1f} %\n"
            f"  η_mech      : {self.n_mech*100:.1f} %\n"
            f"  Work        : {self.work:.2f} W\n"
            f"  Exergy dest.: {self.Ex_D}\n"
            f"  Solved      : {self.Solution_Status}"
        )