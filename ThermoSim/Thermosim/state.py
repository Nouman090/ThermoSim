"""
state.py
--------
The Prop class: represents a single thermodynamic state point.
Given any two independent properties (P, T, H, S, Q, D) and a fluid name
it calculates all the remaining properties via CoolProp.
"""

import warnings
import json
import os
import numpy as np
import CoolProp.CoolProp as CP

from . import config                       # reads dead_states from config.py

# Path to T66.json (sits next to this file OR in project root)
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
_T66_PATHS = [
    os.path.join(_PKG_DIR, 'T66.json'),    # inside the package folder
    os.path.join(_PKG_DIR, '..', 'T66.json'),  # project root
]


class Prop:
    """One thermodynamic state point."""

    _VALID_PROPS = {'P', 'T', 'H', 'S', 'Q', 'D'}

    # ------------------------------------------------------------------ #
    #  constructor
    # ------------------------------------------------------------------ #
    def __init__(self, fluid, StatePointName,
                 Mass_flowrate=None, Solution_Status=False, **properties):
        self.fluid = fluid
        self.Mass_flowrate = Mass_flowrate
        self.StatePointName = StatePointName
        self.Solution_Status = Solution_Status
        self.DeadStates = config.dead_states          # <-- from config

        # --- unwrap the legacy 'pro' wrapper if present -----------------
        if 'pro' in properties:
            properties = properties['pro']

        # keep only non-None values
        self.properties = {k: v for k, v in properties.items()
                           if v is not None}

        # validate keys
        for key in self.properties:
            if key not in self._VALID_PROPS:
                raise ValueError(
                    f"Invalid property '{key}'. "
                    f"Valid keys: {self._VALID_PROPS}"
                )

        if len(self.properties) > 2:
            raise ValueError(
                f"State '{StatePointName}': provide exactly 2 independent "
                f"properties, got {len(self.properties)}: {self.properties}"
            )

        # --- initialise every attribute to None ------------------------
        self.P = self.T = self.H = self.S = None
        self.Q = self.D = self.Cp = None
        self.ex = self.Ex = None

        # set the properties we already know
        for key, value in self.properties.items():
            setattr(self, key, value)

        # --- calculate the rest -----------------------------------------
        if len(self.properties) == 2:
            if self.fluid == 'Therminol66':
                self._calc_therminol66()
            else:
                self._calc_coolprop()

    # ------------------------------------------------------------------ #
    #  CoolProp calculation
    # ------------------------------------------------------------------ #
    def _calc_coolprop(self):
        prop1, prop2 = list(self.properties.keys())
        val1, val2   = list(self.properties.values())

        # --- handle near-saturation for P-T inputs ---------------------
        if {prop1, prop2} == {'P', 'T'}:
            P = self.properties['P']
            T = self.properties['T']
            try:
                T_sat = CP.PropsSI('T', 'P', P, 'Q', 0, self.fluid)
                if abs(T - T_sat) / T_sat < 1e-4:
                    warnings.warn(
                        f"State '{self.StatePointName}' is near saturation. "
                        f"Assuming saturated liquid (Q=0)."
                    )
                    prop1, val1 = 'P', P
                    prop2, val2 = 'Q', 0.0
            except Exception:
                pass

        try:
            if self.P  is None: self.P  = CP.PropsSI('P', prop1, val1, prop2, val2, self.fluid)
            if self.T  is None: self.T  = CP.PropsSI('T', prop1, val1, prop2, val2, self.fluid)
            if self.H  is None: self.H  = CP.PropsSI('H', prop1, val1, prop2, val2, self.fluid)
            if self.S  is None: self.S  = CP.PropsSI('S', prop1, val1, prop2, val2, self.fluid)
            if self.Q  is None: self.Q  = CP.PropsSI('Q', prop1, val1, prop2, val2, self.fluid)
            if self.D  is None: self.D  = CP.PropsSI('D', prop1, val1, prop2, val2, self.fluid)
            if self.Cp is None: self.Cp = CP.PropsSI('C', prop1, val1, prop2, val2, self.fluid)

            # label phase nicely
            self._classify_phase()
            # exergy
            self._calc_exergy()

            self.Solution_Status = True

        except Exception as e:
            raise ValueError(
                f"CoolProp error for state '{self.StatePointName}' "
                f"({self.fluid}): {e}"
            )

    # ------------------------------------------------------------------ #
    def _classify_phase(self):
        """Replace the raw Q number with a readable label for
        superheated / sub-cooled states."""
        try:
            H_liq = CP.PropsSI('H', 'P', self.P, 'Q', 0, self.fluid)
            H_vap = CP.PropsSI('H', 'P', self.P, 'Q', 1, self.fluid)
            T_sat = CP.PropsSI('T', 'P', self.P, 'Q', 0, self.fluid)

            if self.H < H_liq:
                dT = T_sat - self.T
                self.Q = f"Sub-cooled ({dT:.2f} K)"
            elif self.H > H_vap:
                dT = self.T - T_sat
                self.Q = f"Superheated ({dT:.2f} K)"
        except Exception:
            pass  # supercritical, etc.

    # ------------------------------------------------------------------ #
    def _calc_exergy(self):
        """Compute specific flow exergy and total exergy rate."""
        if not self.DeadStates or self.H is None or self.S is None:
            self.ex = self.Ex = None
            return
        T0 = self.DeadStates['T0']
        P0 = self.DeadStates['P0']
        try:
            h0 = CP.PropsSI('H', 'T', T0, 'P', P0, self.fluid)
            s0 = CP.PropsSI('S', 'T', T0, 'P', P0, self.fluid)
        except Exception:
            self.ex = self.Ex = None
            return

        self.ex = (self.H - h0) - T0 * (self.S - s0)
        self.Ex = self.Mass_flowrate * self.ex if self.Mass_flowrate else None

    # ------------------------------------------------------------------ #
    #  Therminol-66 polynomial calculation
    # ------------------------------------------------------------------ #
    def _calc_therminol66(self):
        t66_path = None
        for p in _T66_PATHS:
            if os.path.isfile(p):
                t66_path = p
                break
        if t66_path is None:
            raise FileNotFoundError("T66.json not found.")

        with open(t66_path, 'r') as f:
            prop = json.load(f)

        data   = prop['fluids']['T66']
        h_mean = prop['h_normalization']['h_mean']
        h_std  = prop['h_normalization']['h_std']

        coeffs = {k: data[k]['coeffs'] for k in
                  ['density', 'specific_heat', 'enthalpy',
                   'entropy', 'inverse_enthalpy']}

        if self.T is not None:
            pass                   # T already known
        elif self.H is not None:
            H_norm = (self.H - h_mean) / h_std
            self.T = np.polyval(coeffs['inverse_enthalpy'], H_norm)
        else:
            raise ValueError("Therminol66 requires T or H (plus P).")

        self.D  = np.polyval(coeffs['density'],       self.T)
        self.Cp = np.polyval(coeffs['specific_heat'],  self.T)
        self.H  = np.polyval(coeffs['enthalpy'],       self.T)
        self.S  = np.polyval(coeffs['entropy'],        self.T)

        # exergy
        if self.DeadStates:
            T0 = self.DeadStates['T0']
            h0 = np.polyval(coeffs['enthalpy'], T0)
            s0 = np.polyval(coeffs['entropy'],  T0)
            self.ex = (self.H - h0) - T0 * (self.S - s0)
            self.Ex = (self.Mass_flowrate * self.ex
                       if self.Mass_flowrate else None)

        self.Solution_Status = True

    # ------------------------------------------------------------------ #
    #  public helper (used by Separator, etc.)
    # ------------------------------------------------------------------ #
    def calculate_missing_properties(self):
        if self.fluid == 'Therminol66':
            self._calc_therminol66()
        else:
            self._calc_coolprop()

    # ------------------------------------------------------------------ #
    #  pretty printing
    # ------------------------------------------------------------------ #
    def __repr__(self):
        T_str = f"{self.T - 273.15:.2f} °C" if self.T else "N/A"
        P_str = f"{self.P / 1e5:.2f} bar"   if self.P else "N/A"
        return f"<State '{self.StatePointName}': {self.fluid}, {T_str}, {P_str}>"

    def __str__(self):
        try:
            return (
                f"State: {self.StatePointName}\n"
                f"  Fluid     : {self.fluid}\n"
                f"  P         : {self.P:.2f} Pa  ({self.P/1e5:.2f} bar)\n"
                f"  T         : {self.T - 273.15:.2f} °C  ({self.T:.2f} K)\n"
                f"  H         : {self.H:.2f} J/kg\n"
                f"  S         : {self.S:.4f} J/(kg·K)\n"
                f"  Q         : {self.Q}\n"
                f"  Mass flow : {self.Mass_flowrate} kg/s\n"
                f"  ex        : {self.ex} J/kg\n"
                f"  Ex        : {self.Ex} W\n"
            )
        except Exception:
            return (
                f"State: {self.StatePointName}\n"
                f"  Fluid     : {self.fluid}\n"
                f"  P={self.P}, T={self.T}, H={self.H}, S={self.S}, Q={self.Q}\n"
            )