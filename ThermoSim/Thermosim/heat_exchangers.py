"""
heat_exchangers.py
------------------
HeatExchanger and TES (Thermal Energy Storage) components.

The HeatExchanger logic is kept very close to the original so that
existing models continue to work.  Only bugs are fixed:
  • exergy sign:  (In_hot + In_cold) – (Out_hot + Out_cold)
  • mass-flow-rate resolution uses base-class helper
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve, brentq

from .base_component import Component


# ================================================================== #
#   HEAT EXCHANGER
# ================================================================== #
class HeatExchanger(Component):

    def __init__(self, Model, ID, PPT, HEX_type, HeatAdded,
                 Hot_In_state, Hot_Out_state,
                 Cold_In_state, Cold_Out_state,
                 UA=None, effectiveness=None, Q=None,
                 div_N=200, PPT_graph=False, Calculate=False):

        self.HeatAdded      = HeatAdded
        self.Hot_In_state   = Hot_In_state
        self.Hot_Out_state  = Hot_Out_state
        self.Cold_In_state  = Cold_In_state
        self.Cold_Out_state = Cold_Out_state
        self.HEX_type       = HEX_type
        self.div_N          = div_N
        self.PPT            = PPT
        self.PPT_graph      = PPT_graph
        self.Hot_to_Cold    = None
        self.Q              = Q
        self.UA             = UA
        self.effectiveness  = effectiveness
        self.Hot_Mass_flowrate  = None
        self.Cold_Mass_flowrate = None

        if HeatAdded not in (True, False, None):
            raise ValueError(
                "HeatAdded must be True (adds heat to cycle), "
                "False (rejects heat), or None (internal HEX)."
            )

        # Do NOT call super().__init__ with Calculate yet
        # because we store it and conditionally call Cal()
        self.Model = Model
        self.ID = ID
        self.Solution_Status = False
        self.Ex_D = "Not Calculated"
        self.Model.Component[ID] = self

        if Calculate:
            self.Cal()

    # -------------------------------------------------------------- #
    #  helpers for pinch checks
    # -------------------------------------------------------------- #
    def PPT_Hot_In(self):
        diff = self.Hot_In.T - self.Cold_Out.T
        if abs(diff - self.PPT) < 0.002:
            self.PPT -= 0.001
        elif diff < self.PPT:
            print(f"[{self.ID}] Pinch at hot-inlet side (ΔT = {diff:.2f} K)")

    def PPT_Hot_Out(self):
        diff = self.Hot_Out.T - self.Cold_In.T
        if abs(diff - self.PPT) < 0.002:
            self.PPT -= 0.001
        elif diff < self.PPT:
            raise ValueError(
                f"Pinch at hot-outlet of {self.ID}: ΔT = {diff:.2f} K"
            )

    # -------------------------------------------------------------- #
    #  mass-flow resolution for BOTH sides
    # -------------------------------------------------------------- #
    def _resolve_side(self, pt_in, pt_out, side_name):
        """Return mass-flow rate for one side (hot or cold)."""
        m_in  = pt_in.Mass_flowrate  if pt_in  else None
        m_out = pt_out.Mass_flowrate if pt_out else None

        if m_out is None and m_in is not None:
            if pt_out: pt_out.Mass_flowrate = m_in
            return m_in
        if m_in is None and m_out is not None:
            if pt_in: pt_in.Mass_flowrate = m_out
            return m_out
        if m_in == m_out:
            return m_in
        if m_in is not None and m_out is not None and m_in != m_out:
            raise ValueError(
                f"{side_name}-side mass-flow mismatch in {self.ID}."
            )
        return None

    # -------------------------------------------------------------- #
    #  MAIN CALCULATION
    # -------------------------------------------------------------- #
    def Cal(self):
        # fetch points (may be None for SimpleHEX)
        self.Hot_In   = self.Model.Point[self.Hot_In_state]   if self.Hot_In_state   else None
        self.Hot_Out  = self.Model.Point[self.Hot_Out_state]  if self.Hot_Out_state  else None
        self.Cold_In  = self.Model.Point[self.Cold_In_state]  if self.Cold_In_state  else None
        self.Cold_Out = self.Model.Point[self.Cold_Out_state] if self.Cold_Out_state else None

        # resolve mass-flow rates
        if self.Hot_In is not None or self.Hot_Out is not None:
            self.Hot_Mass_flowrate = self._resolve_side(
                self.Hot_In, self.Hot_Out, 'Hot')
        if self.Cold_In is not None or self.Cold_Out is not None:
            self.Cold_Mass_flowrate = self._resolve_side(
                self.Cold_In, self.Cold_Out, 'Cold')

        # arrays for pinch analysis
        N  = self.div_N
        Th = np.zeros(N + 1)
        Tc = np.zeros(N + 1)
        dT = np.zeros(N + 1)
        h_h = np.zeros(N + 1)
        h_c = np.zeros(N + 1)

        # ---------------------------------------------------------- #
        #  double_pipe / Condenser / Evaporator
        # ---------------------------------------------------------- #
        if self.HEX_type in ('double_pipe', 'Condenser', 'Evaporator'):
            delta_P_c = self.Cold_In.P - self.Cold_Out.P
            delta_P_h = self.Hot_In.P  - self.Hot_Out.P

            if self.Hot_In.H is not None and self.Hot_Out.H is not None:
                if self.Hot_Mass_flowrate is not None and self.Cold_Mass_flowrate is not None:
                    Q = (self.Hot_In.H - self.Hot_Out.H) * self.Hot_Mass_flowrate

                    if self.Cold_Out.H is None:
                        HCO = self.Cold_In.H + Q / self.Cold_Mass_flowrate
                        self.Cold_Out = self.Model.Prop(
                            self.Cold_Out.fluid,
                            StatePointName=self.Cold_Out.StatePointName,
                            H=HCO, P=self.Cold_Out.P,
                            Mass_flowrate=self.Cold_Mass_flowrate)
                        self.PPT_Hot_Out()

                    elif self.Cold_In.H is None:
                        HCI = self.Cold_Out.H - Q / self.Cold_Mass_flowrate
                        self.Cold_In = self.Model.Prop(
                            self.Cold_In.fluid,
                            StatePointName=self.Cold_In.StatePointName,
                            H=HCI, P=self.Cold_In.P,
                            Mass_flowrate=self.Cold_Mass_flowrate)
                        self.PPT_Hot_In()

                    self._fill_pinch_arrays(Th, Tc, dT, h_h, h_c,
                                            Q, delta_P_h, delta_P_c,
                                            direction='cold_to_hot')
                else:
                    # mass-flow unknown → solve via PPT
                    if self.Cold_Out.H is None:
                        self.PPT_Hot_Out()
                        self._solve_unknown_cold_out(
                            Th, Tc, dT, h_h, h_c, delta_P_h, delta_P_c)
                    elif self.Cold_In.H is None:
                        self.PPT_Hot_In()
                        self._solve_unknown_cold_in(
                            Th, Tc, dT, h_h, h_c, delta_P_h, delta_P_c)

            elif self.Cold_In.H is not None and self.Cold_Out.H is not None:
                if self.Hot_Mass_flowrate is not None and self.Cold_Mass_flowrate is not None:
                    Q = (self.Cold_Out.H - self.Cold_In.H) * self.Cold_Mass_flowrate
                    if self.Hot_In.H is not None:
                        HHO = self.Hot_In.H - Q / self.Hot_Mass_flowrate
                        self.Hot_Out = self.Model.Prop(
                            self.Hot_In.fluid,
                            StatePointName=self.Hot_Out.StatePointName,
                            H=HHO, P=self.Hot_Out.P,
                            Mass_flowrate=self.Hot_Mass_flowrate)
                        self.PPT_Hot_In()
                    elif self.Hot_Out.H is not None:
                        HHI = self.Hot_Out.H + Q / self.Hot_Mass_flowrate
                        self.Hot_In = self.Model.Prop(
                            self.Hot_In.fluid,
                            StatePointName=self.Hot_In.StatePointName,
                            H=HHI, P=self.Hot_In.P,
                            Mass_flowrate=self.Hot_Mass_flowrate)
                        self.PPT_Hot_Out()

                    self._fill_pinch_arrays(Th, Tc, dT, h_h, h_c,
                                            Q, delta_P_h, delta_P_c,
                                            direction='hot_to_cold')
                else:
                    if self.Hot_Out.H is None:
                        self.PPT_Hot_In()
                        self._solve_unknown_hot_out(
                            Th, Tc, dT, h_h, h_c, delta_P_h, delta_P_c)
                    elif self.Hot_In.H is None:
                        self.PPT_Hot_Out()
                        self._solve_unknown_hot_in(
                            Th, Tc, dT, h_h, h_c, delta_P_h, delta_P_c)

            elif (self.Cold_In.H is not None and self.Hot_In.H is not None
                  and self.effectiveness is not None
                  and self.Cold_Out.H is None and self.Hot_Out.H is None):
                if self.Hot_Mass_flowrate and self.Cold_Mass_flowrate:
                    self._solve_effectiveness()
                else:
                    raise ValueError(
                        f"Both mass-flow rates needed for effectiveness "
                        f"method in {self.ID}."
                    )

        # ---------------------------------------------------------- #
        #  SimpleHEX  (one-side only)
        # ---------------------------------------------------------- #
        elif self.HEX_type == 'SimpleHEX':
            self._solve_simple_hex()

        else:
            raise ValueError(
                f"Invalid HEX_type '{self.HEX_type}'. "
                f"Use: Evaporator, Condenser, double_pipe, SimpleHEX"
            )

        # ---------------------------------------------------------- #
        #  post-processing (mass-flow ratio, Q, UA, exergy)
        # ---------------------------------------------------------- #
        if self.HEX_type != 'SimpleHEX':
            self.Hot_to_Cold = ((self.Cold_Out.H - self.Cold_In.H)
                                / (self.Hot_In.H  - self.Hot_Out.H))

            if self.Cold_Mass_flowrate is None and self.Hot_Mass_flowrate is not None:
                self.Cold_Mass_flowrate = self.Hot_Mass_flowrate / self.Hot_to_Cold
                self.Cold_In.Mass_flowrate  = self.Cold_Mass_flowrate
                self.Cold_Out.Mass_flowrate = self.Cold_Mass_flowrate
            elif self.Hot_Mass_flowrate is None and self.Cold_Mass_flowrate is not None:
                self.Hot_Mass_flowrate = self.Cold_Mass_flowrate * self.Hot_to_Cold
                self.Hot_In.Mass_flowrate  = self.Hot_Mass_flowrate
                self.Hot_Out.Mass_flowrate = self.Hot_Mass_flowrate

        # Q
        if self.Hot_Mass_flowrate is not None and self.Hot_In is not None and self.Hot_Out is not None:
            self.Q = (self.Hot_In.H - self.Hot_Out.H) * self.Hot_Mass_flowrate
        elif self.Cold_Mass_flowrate is not None and self.Cold_In is not None and self.Cold_Out is not None:
            self.Q = (self.Cold_Out.H - self.Cold_In.H) * self.Cold_Mass_flowrate

        # UA
        if self.HEX_type != 'SimpleHEX' and self.Hot_In and self.Cold_Out and self.Hot_Out and self.Cold_In:
            dT1 = self.Hot_In.T  - self.Cold_Out.T
            dT2 = self.Hot_Out.T - self.Cold_In.T
            if dT1 > 0 and dT2 > 0:
                LMTD = (dT1 - dT2) / np.log(dT1 / dT2) if dT1 != dT2 else dT1
                self.UA = self.Q / LMTD

        # Exergy  (FIXED sign)
        if self.HEX_type != 'SimpleHEX':
            try:
                self.Ex_D = ((self.Hot_In.Ex + self.Cold_In.Ex)
                             - (self.Hot_Out.Ex + self.Cold_Out.Ex))
            except Exception:
                self.Ex_D = "Not Calculated"
        else:
            self.Ex_D = "SimpleHEX: Ex_D not calculated"

        self.Solution_Status = True

        # write back
        if self.Hot_In_state:   self.Model.Point[self.Hot_In_state]   = self.Hot_In
        if self.Hot_Out_state:  self.Model.Point[self.Hot_Out_state]  = self.Hot_Out
        if self.Cold_In_state:  self.Model.Point[self.Cold_In_state]  = self.Cold_In
        if self.Cold_Out_state: self.Model.Point[self.Cold_Out_state] = self.Cold_Out

    # ============================================================== #
    #  internal helper: fill pinch arrays (common pattern)
    # ============================================================== #
    def _fill_pinch_arrays(self, Th, Tc, dT, h_h, h_c,
                           Q, dP_h, dP_c, direction='cold_to_hot'):
        N = self.div_N
        q = Q / N

        if direction == 'cold_to_hot':
            Th[0]  = self.Hot_Out.T;  h_h[0] = self.Hot_Out.H
            Tc[0]  = self.Cold_In.T;  h_c[0] = self.Cold_In.H
            for n in range(1, N + 1):
                h_h[n] = h_h[n-1] + q / self.Hot_Mass_flowrate
                h_c[n] = h_c[n-1] + q / self.Cold_Mass_flowrate
                Ph = self.Hot_Out.P  + (dP_h / N) * n
                Pc = self.Cold_In.P  - (dP_c / N) * n
                Th[n] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h', H=h_h[n], P=Ph).T
                Tc[n] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c', H=h_c[n], P=Pc).T
                dT[n] = Th[n] - Tc[n]
            dT[0] = Th[0] - Tc[0]

        else:  # hot_to_cold
            Th[0]  = self.Hot_In.T;   h_h[0] = self.Hot_In.H
            Tc[0]  = self.Cold_Out.T; h_c[0] = self.Cold_Out.H
            for n in range(1, N + 1):
                h_h[n] = h_h[n-1] - q / self.Hot_Mass_flowrate
                h_c[n] = h_c[n-1] - q / self.Cold_Mass_flowrate
                Ph = self.Hot_In.P  - (dP_h / N) * n
                Pc = self.Cold_Out.P + (dP_c / N) * n
                Th[n] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h', H=h_h[n], P=Ph).T
                Tc[n] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c', H=h_c[n], P=Pc).T
                dT[n] = Th[n] - Tc[n]
            dT[0] = Th[0] - Tc[0]

        if self.PPT_graph:
            self._draw_pinch(Th, Tc, dT)

    # ============================================================== #
    #  PPT-based solvers (brentq) — kept close to original logic
    # ============================================================== #
    def _solve_unknown_cold_out(self, Th, Tc, dT, h_h, h_c, dP_h, dP_c):
        N = self.div_N

        def f(T):
            Q_sp = self.Hot_In.H - self.Hot_Out.H
            if self.HEX_type in ('double_pipe', 'Condenser'):
                self.Cold_Out = self.Model.Prop(
                    self.Cold_Out.fluid,
                    StatePointName=self.Cold_Out.StatePointName,
                    T=T, P=self.Cold_Out.P)
            elif self.HEX_type == 'Evaporator':
                T_sat = self.Model.Prop(
                    self.Cold_In.fluid, StatePointName='_sat',
                    P=self.Cold_Out.P, Q=1).T
                self.Cold_Out = self.Model.Prop(
                    self.Cold_In.fluid,
                    StatePointName=self.Cold_Out.StatePointName,
                    T=T_sat + T, P=self.Cold_Out.P)

            m_c = Q_sp / (self.Cold_Out.H - self.Cold_In.H)
            q = Q_sp / N
            Th[0] = self.Hot_Out.T;  h_h[0] = self.Hot_Out.H
            Tc[0] = self.Cold_In.T;  h_c[0] = self.Cold_In.H
            dT[0] = Th[0] - Tc[0]
            for n in range(1, N + 1):
                h_h[n] = h_h[n-1] + q
                h_c[n] = h_c[n-1] + q / m_c
                Ph = self.Hot_Out.P + (dP_h / N) * n
                Pc = self.Cold_In.P - (dP_c / N) * n
                Th[n] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h', H=h_h[n], P=Ph).T
                Tc[n] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c', H=h_c[n], P=Pc).T
                dT[n] = Th[n] - Tc[n]
            return min(dT) - self.PPT

        try:
            if self.HEX_type in ('double_pipe', 'Condenser'):
                T_co = brentq(f, self.Cold_In.T + 1e-7, self.Hot_In.T,
                              xtol=1e-3, rtol=1e-3)
                self.Cold_Out = self.Model.Prop(
                    self.Cold_In.fluid,
                    StatePointName=self.Cold_Out.StatePointName,
                    T=T_co, P=self.Cold_Out.P,
                    Mass_flowrate=self.Cold_Mass_flowrate)
            elif self.HEX_type == 'Evaporator':
                T_sat = self.Model.Prop(
                    self.Cold_Out.fluid, StatePointName='_sat',
                    Q=1, P=self.Cold_Out.P).T
                T_sup = brentq(f, 0.001,
                               CP.PropsSI('TMAX', self.Cold_Out.fluid) - T_sat,
                               xtol=1e-3, rtol=1e-3)
                self.Cold_Out = self.Model.Prop(
                    self.Cold_In.fluid,
                    StatePointName=self.Cold_Out.StatePointName,
                    T=T_sat + T_sup, P=self.Cold_Out.P)

            if self.PPT_graph:
                self._draw_pinch(Th, Tc, dT)
        except Exception as e:
            print(f"[{self.ID}] solver error: {e}")
            if self.PPT_graph:
                self._draw_pinch(Th, Tc, dT)

    def _solve_unknown_cold_in(self, Th, Tc, dT, h_h, h_c, dP_h, dP_c):
        N = self.div_N

        def f(T):
            Q_sp = self.Hot_In.H - self.Hot_Out.H
            self.Cold_In = self.Model.Prop(
                self.Cold_In.fluid,
                StatePointName=self.Cold_In.StatePointName,
                T=T, P=self.Cold_In.P)
            m_c = Q_sp / (self.Cold_Out.H - self.Cold_In.H)
            q = Q_sp / N
            Th[0] = self.Hot_Out.T;  h_h[0] = self.Hot_Out.H
            Tc[0] = self.Cold_In.T;  h_c[0] = self.Cold_In.H
            dT[0] = Th[0] - Tc[0]
            for n in range(1, N + 1):
                h_h[n] = h_h[n-1] + q
                h_c[n] = h_c[n-1] + q / m_c
                Ph = self.Hot_Out.P + (dP_h / N) * n
                Pc = self.Cold_In.P - (dP_c / N) * n
                Th[n] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h', H=h_h[n], P=Ph).T
                Tc[n] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c', H=h_c[n], P=Pc).T
                dT[n] = Th[n] - Tc[n]
            return min(dT) - self.PPT

        try:
            T_ci = brentq(f,
                          CP.PropsSI('TMIN', self.Cold_In.fluid),
                          self.Hot_In.T,
                          xtol=1e-3, rtol=1e-3)
            self.Cold_In = self.Model.Prop(
                self.Cold_In.fluid,
                StatePointName=self.Cold_In.StatePointName,
                T=T_ci, P=self.Cold_In.P,
                Mass_flowrate=self.Cold_Mass_flowrate)
            if self.PPT_graph:
                self._draw_pinch(Th, Tc, dT)
        except Exception as e:
            print(f"[{self.ID}] solver error: {e}")

    def _solve_unknown_hot_out(self, Th, Tc, dT, h_h, h_c, dP_h, dP_c):
        N = self.div_N

        def f(T):
            Q_sp = self.Cold_Out.H - self.Cold_In.H
            if self.HEX_type in ('double_pipe', 'Evaporator'):
                self.Hot_Out = self.Model.Prop(
                    self.Hot_In.fluid,
                    StatePointName=self.Hot_Out.StatePointName,
                    T=T, P=self.Hot_Out.P,
                    Mass_flowrate=self.Hot_Mass_flowrate)
            elif self.HEX_type == 'Condenser':
                T_sat = self.Model.Prop(
                    self.Hot_In.fluid, StatePointName='_sat',
                    P=self.Hot_Out.P, Q=1).T
                self.Hot_Out = self.Model.Prop(
                    self.Hot_In.fluid,
                    StatePointName=self.Hot_Out.StatePointName,
                    T=T_sat - T, P=self.Hot_Out.P,
                    Mass_flowrate=self.Hot_Mass_flowrate)

            m_h = Q_sp / (self.Hot_In.H - self.Hot_Out.H)
            q = Q_sp / N
            Th[0] = self.Hot_In.T;   h_h[0] = self.Hot_In.H
            Tc[0] = self.Cold_Out.T; h_c[0] = self.Cold_Out.H
            dT[0] = Th[0] - Tc[0]
            for n in range(1, N + 1):
                h_h[n] = h_h[n-1] - q / m_h
                h_c[n] = h_c[n-1] - q
                Ph = self.Hot_In.P  - (dP_h / N) * n
                Pc = self.Cold_Out.P + (dP_c / N) * n
                Th[n] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h', H=h_h[n], P=Ph).T
                Tc[n] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c', H=h_c[n], P=Pc).T
                dT[n] = Th[n] - Tc[n]
            return min(dT) - self.PPT

        try:
            if self.HEX_type in ('double_pipe', 'Evaporator'):
                T_ho = brentq(f, self.Cold_In.T, self.Hot_In.T - 0.001,
                              xtol=1e-3, rtol=1e-3)
                self.Hot_Out = self.Model.Prop(
                    self.Hot_In.fluid,
                    StatePointName=self.Hot_Out.StatePointName,
                    T=T_ho, P=self.Hot_Out.P,
                    Mass_flowrate=self.Hot_Mass_flowrate)
            elif self.HEX_type == 'Condenser':
                T_sup = brentq(f, 0.1, 100, xtol=1e-3, rtol=1e-3)
                T_ho = self.Model.Prop(
                    self.Hot_Out.fluid, StatePointName='_sat',
                    Q=1, P=self.Hot_Out.P).T - T_sup
                self.Hot_Out = self.Model.Prop(
                    self.Hot_Out.fluid,
                    StatePointName=self.Hot_Out.StatePointName,
                    T=T_ho, P=self.Hot_Out.P,
                    Mass_flowrate=self.Hot_Mass_flowrate)
            if self.PPT_graph:
                self._draw_pinch(Th, Tc, dT)
        except Exception as e:
            print(f"[{self.ID}] solver error: {e}")

    def _solve_unknown_hot_in(self, Th, Tc, dT, h_h, h_c, dP_h, dP_c):
        N = self.div_N

        def f(T):
            Q_sp = self.Cold_Out.H - self.Cold_In.H
            self.Hot_In = self.Model.Prop(
                self.Hot_In.fluid,
                StatePointName=self.Hot_In.StatePointName,
                T=T, P=self.Hot_In.P,
                Mass_flowrate=self.Hot_Mass_flowrate)
            m_h = Q_sp / (self.Hot_In.H - self.Hot_Out.H)
            q = Q_sp / N
            Th[0] = self.Hot_In.T;   h_h[0] = self.Hot_In.H
            Tc[0] = self.Cold_Out.T; h_c[0] = self.Cold_Out.H
            dT[0] = Th[0] - Tc[0]
            for n in range(1, N + 1):
                h_h[n] = h_h[n-1] - q / m_h
                h_c[n] = h_c[n-1] - q
                Ph = self.Hot_In.P   - (dP_h / N) * n
                Pc = self.Cold_Out.P + (dP_c / N) * n
                Th[n] = self.Model.Prop(self.Hot_In.fluid,  StatePointName='_h', H=h_h[n], P=Ph).T
                Tc[n] = self.Model.Prop(self.Cold_In.fluid, StatePointName='_c', H=h_c[n], P=Pc).T
                dT[n] = Th[n] - Tc[n]
            return min(dT) - self.PPT

        try:
            T_hi = brentq(f,
                          self.Hot_Out.T + 0.001,
                          CP.PropsSI('TMAX', self.Hot_In.fluid),
                          xtol=1e-3, rtol=1e-3)
            self.Hot_In = self.Model.Prop(
                self.Hot_In.fluid,
                StatePointName=self.Hot_In.StatePointName,
                T=T_hi, P=self.Hot_In.P,
                Mass_flowrate=self.Hot_Mass_flowrate)
            if self.PPT_graph:
                self._draw_pinch(Th, Tc, dT)
        except Exception as e:
            print(f"[{self.ID}] solver error: {e}")

    # ============================================================== #
    #  effectiveness method
    # ============================================================== #
    def _solve_effectiveness(self):
        Ch = self.Hot_Mass_flowrate  * self.Hot_In.Cp
        Cc = self.Cold_Mass_flowrate * self.Cold_In.Cp
        for _ in range(50):
            Cmin = min(Ch, Cc)
            self.Q = self.effectiveness * Cmin * (self.Hot_In.T - self.Cold_In.T)
            hh = self.Hot_In.H  - self.Q / self.Hot_Mass_flowrate
            hc = self.Cold_In.H + self.Q / self.Cold_Mass_flowrate
            self.Hot_Out = self.Model.Prop(
                self.Hot_Out.fluid,
                StatePointName=self.Hot_Out.StatePointName,
                P=self.Hot_Out.P, H=hh,
                Mass_flowrate=self.Hot_Mass_flowrate)
            self.Cold_Out = self.Model.Prop(
                self.Cold_Out.fluid,
                StatePointName=self.Cold_Out.StatePointName,
                P=self.Cold_Out.P, H=hc,
                Mass_flowrate=self.Cold_Mass_flowrate)
            ch_avg = (self.Hot_In.Cp  + self.Hot_Out.Cp)  * self.Hot_Mass_flowrate  / 2
            cc_avg = (self.Cold_In.Cp + self.Cold_Out.Cp) * self.Cold_Mass_flowrate / 2
            if abs(Ch - ch_avg) / ch_avg < 1e-3 and abs(Cc - cc_avg) / cc_avg < 1e-3:
                break
            Ch, Cc = ch_avg, cc_avg

    # ============================================================== #
    #  SimpleHEX solver
    # ============================================================== #
    def _solve_simple_hex(self):
        CI = self.Cold_In
        CO = self.Cold_Out
        HI = self.Hot_In
        HO = self.Hot_Out
        mc = self.Cold_Mass_flowrate
        mh = self.Hot_Mass_flowrate

        if CI and CO is not None and CO.H is None and self.Q and mc:
            H = self.Q / mc + CI.H
            self.Cold_Out = self.Model.Prop(
                CO.fluid, StatePointName=CO.StatePointName,
                P=CO.P, H=H, Mass_flowrate=mc)

        elif CO and CI is not None and CI.H is None and self.Q and mc:
            H = CO.H - self.Q / mc
            self.Cold_In = self.Model.Prop(
                CI.fluid, StatePointName=CI.StatePointName,
                P=CI.P, H=H, Mass_flowrate=mc)

        elif CI and CO and self.Q is None and mc:
            self.Q = (CO.H - CI.H) * mc

        elif HI and HO is not None and HO.H is None and self.Q and mh:
            H = HI.H - self.Q / mh
            self.Hot_Out = self.Model.Prop(
                HO.fluid, StatePointName=HO.StatePointName,
                P=HO.P, H=H, Mass_flowrate=mh)

        elif HO and HI is not None and HI.H is None and self.Q and mh:
            H = HO.H + self.Q / mh
            self.Hot_In = self.Model.Prop(
                HI.fluid, StatePointName=HI.StatePointName,
                P=HI.P, H=H, Mass_flowrate=mh)

        elif HI and HO and self.Q is None and mh:
            self.Q = (HI.H - HO.H) * mh

        else:
            raise ValueError(
                f"SimpleHEX {self.ID}: need mass-flow rate and 2 of 3 "
                f"(inlet, outlet, Q) on one side."
            )

    # ============================================================== #
    #  pinch-point plot
    # ============================================================== #
    def _draw_pinch(self, Th, Tc, dT):
        plt.figure()
        plt.plot(range(len(Th)), Th - 273.15, 'r-', label='Hot')
        plt.plot(range(len(Tc)), Tc - 273.15, 'b-', label='Cold')
        plt.legend()
        idx = np.argmin(dT)
        plt.title(f"{self.ID}  —  min ΔT = {min(dT):.2f} K")
        plt.xlabel('Segment')
        plt.ylabel('Temperature [°C]')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    # ============================================================== #
    def __str__(self):
        if self.HEX_type == 'SimpleHEX':
            if self.Cold_In and self.Hot_In is None:
                T_in  = self.Cold_In.T  - 273.15
                T_out = self.Cold_Out.T - 273.15
            elif self.Hot_In:
                T_in  = self.Hot_In.T  - 273.15
                T_out = self.Hot_Out.T - 273.15
            else:
                T_in = T_out = '?'
            return (
                f"{self.ID} (SimpleHEX):\n"
                f"  T_in  : {T_in} °C\n"
                f"  T_out : {T_out} °C\n"
                f"  Q     : {self.Q} W\n"
                f"  Solved: {self.Solution_Status}"
            )
        try:
            return (
                f"{self.ID} (HEX - {self.HEX_type}):\n"
                f"  Hot  : {self.Hot_In.T-273.15:.1f} → {self.Hot_Out.T-273.15:.1f} °C"
                f"  ({self.Hot_Mass_flowrate} kg/s)\n"
                f"  Cold : {self.Cold_In.T-273.15:.1f} → {self.Cold_Out.T-273.15:.1f} °C"
                f"  ({self.Cold_Mass_flowrate} kg/s)\n"
                f"  Q    : {self.Q} W\n"
                f"  UA   : {self.UA} W/K\n"
                f"  Ex_D : {self.Ex_D}\n"
                f"  Solved: {self.Solution_Status}"
            )
        except Exception:
            return f"{self.ID}: solved={self.Solution_Status}"


# ================================================================== #
#   TES  (Thermal Energy Storage)
# ================================================================== #
class TES(Component):

    def __init__(self, Model, ID, PPT, Charge, T_melt,
                 Hot_In_state, Hot_Out_state,
                 Cold_In_state, Cold_Out_state,
                 Charging_time, Discharging_time,
                 per_loss, Capacity=None, Calculate=False):
        self.T_melt           = T_melt
        self.per_loss         = per_loss
        self.Hot_In_state     = Hot_In_state
        self.Hot_Out_state    = Hot_Out_state
        self.Cold_In_state    = Cold_In_state
        self.Cold_Out_state   = Cold_Out_state
        self.Charging_time    = Charging_time
        self.Discharging_time = Discharging_time
        self.Charge           = Charge
        self.Capacity         = Capacity
        self.PPT              = PPT
        self.Charging_Power   = None
        self.Discharging_Power = None

        # manual registration (don't call super().__init__ with Calculate
        # because we need all attributes set first)
        self.Model = Model
        self.ID = ID
        self.Solution_Status = False
        self.Ex_D = "Not Calculated"
        self.Model.Component[ID] = self

        if Calculate:
            self.Cal()

    def Cal(self):
        self.Hot_In   = self.Model.Point[self.Hot_In_state]   if self.Hot_In_state   else None
        self.Hot_Out  = self.Model.Point[self.Hot_Out_state]  if self.Hot_Out_state  else None
        self.Cold_In  = self.Model.Point[self.Cold_In_state]  if self.Cold_In_state  else None
        self.Cold_Out = self.Model.Point[self.Cold_Out_state] if self.Cold_Out_state else None

        # --- resolve mass flows ---
        self.Hot_Mass_flowrate = None
        self.Cold_Mass_flowrate = None

        if self.Hot_In is not None and self.Hot_Out is not None:
            mi = self.Hot_In.Mass_flowrate
            mo = self.Hot_Out.Mass_flowrate
            if mo is None and mi is not None:
                self.Hot_Out.Mass_flowrate = mi; self.Hot_Mass_flowrate = mi
            elif mi is None and mo is not None:
                self.Hot_In.Mass_flowrate = mo;  self.Hot_Mass_flowrate = mo
            elif mi == mo:
                self.Hot_Mass_flowrate = mi
            else:
                raise ValueError(f"Hot mass-flow mismatch in {self.ID}")

        if self.Cold_In is not None and self.Cold_Out is not None:
            mi = self.Cold_In.Mass_flowrate
            mo = self.Cold_Out.Mass_flowrate
            if mo is None and mi is not None:
                self.Cold_Out.Mass_flowrate = mi; self.Cold_Mass_flowrate = mi
            elif mi is None and mo is not None:
                self.Cold_In.Mass_flowrate = mo;  self.Cold_Mass_flowrate = mo
            elif mi == mo:
                self.Cold_Mass_flowrate = mi
            else:
                raise ValueError(f"Cold mass-flow mismatch in {self.ID}")

        # --- Discharging ---
        if self.Charge == 'Discharging':
            if self.Capacity is not None:
                self.CapacityD = self.Capacity * (1 - self.per_loss)
            else:
                self.CapacityD = None

            CI = self.Cold_In
            CO = self.Cold_Out
            mc = self.Cold_Mass_flowrate
            dt = self.Discharging_time * 3600

            if CI.H is not None and CO.H is None and self.CapacityD is None:
                if self.PPT < (self.T_melt - CI.T):
                    self.Cold_Out = self.Model.Prop(
                        CO.fluid, StatePointName=CO.StatePointName,
                        P=CO.P, T=self.T_melt - self.PPT,
                        Mass_flowrate=mc)
                    self.CapacityD = (self.Cold_Out.H - CI.H) * mc * dt
                else:
                    raise ValueError(f"(T_melt - PPT) < T_in in {self.ID}")

            elif CI.H is not None and CO.H is None and self.CapacityD is not None:
                H = CI.H + self.CapacityD / (mc * dt)
                self.Cold_Out = self.Model.Prop(
                    CO.fluid, StatePointName=CO.StatePointName,
                    P=CO.P, H=H, Mass_flowrate=mc)

            elif CI.H is not None and CO.H is not None and self.CapacityD is None:
                self.CapacityD = (CO.H - CI.H) * mc * dt

            elif CI.H is None and CO.H is None and self.CapacityD is not None:
                self.Cold_Out = self.Model.Prop(
                    CO.fluid, StatePointName=CO.StatePointName,
                    P=CO.P, T=self.T_melt - self.PPT, Mass_flowrate=mc)
                H = self.Cold_Out.H - self.CapacityD / (mc * dt)
                self.Cold_In = self.Model.Prop(
                    CI.fluid, StatePointName=CI.StatePointName,
                    P=CI.P, H=H, Mass_flowrate=mc)

            elif CI.H is None and CO.H is not None and self.CapacityD is not None:
                H = CO.H - self.CapacityD / (mc * dt)
                self.Cold_In = self.Model.Prop(
                    CI.fluid, StatePointName=CI.StatePointName,
                    P=CI.P, H=H, Mass_flowrate=mc)

            self.Discharging_Power = self.CapacityD / dt

        # --- Charging ---
        elif self.Charge == 'Charging':
            HI = self.Hot_In
            HO = self.Hot_Out
            mh = self.Hot_Mass_flowrate
            dt = self.Charging_time * 3600

            if HI.H is not None and HO.H is None and self.Capacity is None:
                if self.PPT < (HI.T - self.T_melt):
                    self.Hot_Out = self.Model.Prop(
                        HO.fluid, StatePointName=HO.StatePointName,
                        P=HO.P, T=self.T_melt + self.PPT,
                        Mass_flowrate=mh)
                    self.Capacity = (HI.H - self.Hot_Out.H) * mh * dt
                else:
                    raise ValueError(f"(T_melt + PPT) > T_in in {self.ID}")

            elif HI.H is not None and HO.H is None and self.Capacity is not None:
                H = HI.H - self.Capacity / (mh * dt)
                self.Hot_Out = self.Model.Prop(
                    HO.fluid, StatePointName=HO.StatePointName,
                    P=HO.P, H=H, Mass_flowrate=mh)

            elif HI.H is not None and HO.H is not None and self.Capacity is None:
                self.Capacity = (HI.H - HO.H) * mh * dt

            elif HI.H is None and HO.H is None and self.Capacity is not None:
                self.Hot_Out = self.Model.Prop(
                    HO.fluid, StatePointName=HO.StatePointName,
                    P=HO.P, T=self.T_melt + self.PPT, Mass_flowrate=mh)
                H = self.Hot_Out.H + self.Capacity / (mh * dt)
                self.Hot_In = self.Model.Prop(
                    HI.fluid, StatePointName=HI.StatePointName,
                    P=HI.P, H=H, Mass_flowrate=mh)

            elif HI.H is None and HO.H is not None and self.Capacity is not None:
                H = HO.H + self.Capacity / (mh * dt)
                self.Hot_In = self.Model.Prop(
                    HI.fluid, StatePointName=HI.StatePointName,
                    P=HI.P, H=H, Mass_flowrate=mh)

            self.Charging_Power = self.Capacity / dt
        else:
            raise ValueError(
                f"Invalid Charge='{self.Charge}' in {self.ID}. "
                f"Use 'Charging' or 'Discharging'."
            )

        # exergy (FIXED sign)
        try:
            self.Ex_D = ((self.Hot_In.Ex + self.Cold_In.Ex)
                         - (self.Hot_Out.Ex + self.Cold_Out.Ex))
        except Exception:
            self.Ex_D = "Not Calculated"

        self.Solution_Status = True

        if self.Hot_In_state:   self.Model.Point[self.Hot_In_state]   = self.Hot_In
        if self.Hot_Out_state:  self.Model.Point[self.Hot_Out_state]  = self.Hot_Out
        if self.Cold_In_state:  self.Model.Point[self.Cold_In_state]  = self.Cold_In
        if self.Cold_Out_state: self.Model.Point[self.Cold_Out_state] = self.Cold_Out

    def __str__(self):
        return (
            f"{self.ID} (TES - {self.Charge}):\n"
            f"  Charging Power   : {self.Charging_Power} W\n"
            f"  Discharging Power: {self.Discharging_Power} W\n"
            f"  Capacity         : {self.Capacity} J\n"
            f"  Solved           : {self.Solution_Status}"
        )