"""
state.py
--------
The Prop class: represents a single thermodynamic state point.
Given any two independent properties (P, T, H, S, Q, D) and a fluid name
it calculates all the remaining properties via CoolProp.

Auto-recalculation behaviour
-----------------------------
After construction is complete, assigning to any watched attribute triggers
an automatic recalculation:

  Thermodynamic inputs  → resets ``self.properties`` to the two new values,
                          then reruns the full solve.
  ``Mass_flowrate``     → reruns only ``_calc_exergy`` (cheap update of Ex).
  ``DeadStates``        → same as Mass_flowrate (cheap exergy-only update).

To change TWO inputs atomically (without a redundant intermediate solve) use
the ``update(**kwargs)`` helper:

    S1.update(P=5e6, T=700)
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

# Attributes whose assignment triggers a full thermodynamic re-solve.
_THERMO_INPUTS = {'P', 'T', 'H', 'S', 'Q', 'D'}

# Fluid names handled by the built-in Therminol-66 correlations rather than
# by CoolProp.  Matching is case-insensitive.
_T66_ALIASES = {'therminol66', 'therminol-66', 't66', 'incomp::t66'}

# Attributes whose assignment triggers only the cheap exergy update.
_EXERGY_INPUTS = {'Mass_flowrate', 'DeadStates'}

# Sentinel so we can detect "no propagator registered yet".
_NO_PROPAGATOR = object()


class Prop:
    """One thermodynamic state point."""

    _VALID_PROPS = {'P', 'T', 'H', 'S', 'Q', 'D'}

    # ------------------------------------------------------------------ #
    #  constructor
    # ------------------------------------------------------------------ #
    def __init__(self, fluid, StatePointName,
                 Mass_flowrate=None, Solution_Status=False, **properties):

        # ---- Guard flag: disable __setattr__ magic during __init__ ----
        # We write directly to __dict__ to avoid triggering our own hook.
        object.__setattr__(self, '_initialising', True)

        # Flow-propagation callback.  Set by ThermodynamicModel after the
        # point is registered.  Signature: propagator(point_name, new_value)
        object.__setattr__(self, '_flow_propagator', _NO_PROPAGATOR)

        # Flag: was Mass_flowrate explicitly supplied by the user at
        # construction time?  Used by _replay_known_flowrates() to find
        # the correct seed points even inside closed loops, where every
        # point technically has an upstream neighbour.
        object.__setattr__(self, '_user_set_flowrate', Mass_flowrate is not None)

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
        self.phase = None

        # set the properties we already know
        for key, value in self.properties.items():
            setattr(self, key, value)

        # --- calculate the rest -----------------------------------------
        if len(self.properties) == 2:
            if self._is_therminol66(self.fluid):
                self._calc_therminol66()
            else:
                self._calc_coolprop()

        # ---- Construction done: enable the auto-recalculation hook ----
        object.__setattr__(self, '_initialising', False)

    # ------------------------------------------------------------------ #
    #  auto-recalculation hook
    # ------------------------------------------------------------------ #
    def __setattr__(self, name, value):
        # During __init__ (or internal calc methods) just set normally.
        if self.__dict__.get('_initialising', True):
            object.__setattr__(self, name, value)
            return

        # Always set the value first so the calc methods see the new data.
        object.__setattr__(self, name, value)

        if name in _THERMO_INPUTS:
            # One of the two defining inputs changed.
            # Rebuild self.properties from the current numeric values of
            # P/T/H/S/Q/D.  (The isinstance(str) filter is defensive only —
            # since v3.2.2 Q is always numeric or None.)
            new_props = {}
            for prop in self._VALID_PROPS:
                v = self.__dict__.get(prop)
                if v is not None and not isinstance(v, str):
                    new_props[prop] = v

            if len(new_props) >= 2:
                # The attribute that was just assigned always wins; the
                # second input is then taken from the remaining known values
                # in the fixed precedence order P → T → H → S → Q → D.
                # (This is a precedence rule, not a recency rule: assigning P
                # to a point that knows both T and H re-solves at constant T.)
                ordered = [name] + [k for k in ('P', 'T', 'H', 'S', 'Q', 'D')
                                     if k != name]
                two = {}
                for k in ordered:
                    v = new_props.get(k)
                    if v is not None:
                        two[k] = v
                    if len(two) == 2:
                        break
                object.__setattr__(self, 'properties', two)
            else:
                # Fewer than two inputs known: keep whatever we do have.
                # Without this the reset loop below wipes the value that was
                # just assigned (self.properties is still {}), so
                # ``pt.P = 1e6`` on a fresh point silently left P as None.
                object.__setattr__(self, 'properties', new_props)

            # Reset all derived attributes before re-solving
            for attr in ('P', 'T', 'H', 'S', 'Q', 'D', 'Cp', 'ex', 'Ex',
                         'phase'):
                if attr not in self.properties:
                    object.__setattr__(self, attr, None)
            # Re-apply the defining inputs (they might have been wiped above)
            for k, v in self.properties.items():
                object.__setattr__(self, k, v)

            object.__setattr__(self, 'Solution_Status', False)
            # Only attempt a solve once two independent inputs are present.
            # Assigning the FIRST property to an under-defined point must
            # simply store the value, not raise from CoolProp.
            if len(self.properties) == 2:
                self.calculate_missing_properties()

        elif name in _EXERGY_INPUTS:
            # Only Mass_flowrate or DeadStates changed — re-run exergy only.
            if self.__dict__.get('Solution_Status'):
                self._calc_exergy()

            # If a flow-propagation network has been registered AND this is
            # a Mass_flowrate change (not DeadStates), propagate downstream.
            if name == 'Mass_flowrate':
                propagator = self.__dict__.get('_flow_propagator', _NO_PROPAGATOR)
                if propagator is not _NO_PROPAGATOR and propagator is not None:
                    propagator(self.StatePointName, value)

    # ------------------------------------------------------------------ #
    #  atomic multi-property update (avoids intermediate re-solves)
    # ------------------------------------------------------------------ #
    def update(self, **kwargs):
        """Change one or more defining inputs without triggering
        intermediate recalculations.

        Example
        -------
        S1.update(P=5e6, T=700)   # recalculates once, not twice
        """
        # Suspend the hook during the batch update.
        object.__setattr__(self, '_initialising', True)
        try:
            for k, v in kwargs.items():
                object.__setattr__(self, k, v)
            # Rebuild properties dict from supplied kwargs only (must be 2).
            if not kwargs.keys() <= self._VALID_PROPS | _EXERGY_INPUTS:
                raise ValueError(
                    f"Unknown keys in update(): "
                    f"{kwargs.keys() - self._VALID_PROPS - _EXERGY_INPUTS}"
                )
            thermo_kwargs = {k: v for k, v in kwargs.items()
                             if k in self._VALID_PROPS}
            if thermo_kwargs:
                object.__setattr__(self, 'properties', thermo_kwargs)
                # Reset derived attributes
                for attr in ('P', 'T', 'H', 'S', 'Q', 'D', 'Cp', 'ex', 'Ex',
                             'phase'):
                    if attr not in thermo_kwargs:
                        object.__setattr__(self, attr, None)
                for k, v in thermo_kwargs.items():
                    object.__setattr__(self, k, v)
                object.__setattr__(self, 'Solution_Status', False)
        finally:
            object.__setattr__(self, '_initialising', False)

        if len(self.properties) == 2:
            self.calculate_missing_properties()

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

        # Guard: temporarily disable the hook so internal setattr calls
        # inside this method don't re-trigger recalculation.
        object.__setattr__(self, '_initialising', True)
        try:
            if self.P  is None: object.__setattr__(self, 'P',  CP.PropsSI('P', prop1, val1, prop2, val2, self.fluid))
            if self.T  is None: object.__setattr__(self, 'T',  CP.PropsSI('T', prop1, val1, prop2, val2, self.fluid))
            if self.H  is None: object.__setattr__(self, 'H',  CP.PropsSI('H', prop1, val1, prop2, val2, self.fluid))
            if self.S  is None: object.__setattr__(self, 'S',  CP.PropsSI('S', prop1, val1, prop2, val2, self.fluid))
            if self.Q  is None: object.__setattr__(self, 'Q',  CP.PropsSI('Q', prop1, val1, prop2, val2, self.fluid))
            if self.D  is None: object.__setattr__(self, 'D',  CP.PropsSI('D', prop1, val1, prop2, val2, self.fluid))
            if self.Cp is None: object.__setattr__(self, 'Cp', CP.PropsSI('C', prop1, val1, prop2, val2, self.fluid))

            self._classify_phase()
            self._calc_exergy()

            object.__setattr__(self, 'Solution_Status', True)

        except Exception as e:
            raise ValueError(
                f"CoolProp error for state '{self.StatePointName}' "
                f"({self.fluid}): {e}"
            )
        finally:
            object.__setattr__(self, '_initialising', False)

    # ------------------------------------------------------------------ #
    def _classify_phase(self):
        """
        Set the human-readable ``self.phase`` label.

        ``self.Q`` is ALWAYS left numeric (0-1 inside the dome) or set to
        ``None`` outside it.  CoolProp returns -1 for single-phase states,
        which is not a meaningful quality, so it is normalised to ``None``.

        Historically (<= v3.2.1) this method overwrote ``Q`` with a string
        such as ``"Superheated (46.97 K)"``, which broke every consumer that
        did arithmetic on ``Q`` (Separator, the flow-propagation network,
        save/load).  The label now lives in ``phase`` instead.
        """
        q = self.__dict__.get('Q')

        # Normalise CoolProp's -1 sentinel (and anything outside [0,1]).
        if isinstance(q, (int, float)) and not (0.0 <= q <= 1.0):
            object.__setattr__(self, 'Q', None)
            q = None

        try:
            H_liq = CP.PropsSI('H', 'P', self.P, 'Q', 0, self.fluid)
            H_vap = CP.PropsSI('H', 'P', self.P, 'Q', 1, self.fluid)
            T_sat = CP.PropsSI('T', 'P', self.P, 'Q', 0, self.fluid)

            if self.H < H_liq:
                object.__setattr__(self, 'phase',
                                   f"Sub-cooled ({T_sat - self.T:.2f} K)")
            elif self.H > H_vap:
                object.__setattr__(self, 'phase',
                                   f"Superheated ({self.T - T_sat:.2f} K)")
            elif q == 0.0:
                object.__setattr__(self, 'phase', "Saturated liquid")
            elif q == 1.0:
                object.__setattr__(self, 'phase', "Saturated vapour")
            else:
                object.__setattr__(self, 'phase', f"Two-phase (Q={q:.4f})")
        except Exception:
            object.__setattr__(self, 'phase', "Supercritical / undefined")

    # ------------------------------------------------------------------ #
    def _calc_exergy(self):
        """Compute specific flow exergy and total exergy rate."""
        if not self.DeadStates or self.H is None or self.S is None:
            object.__setattr__(self, 'ex', None)
            object.__setattr__(self, 'Ex', None)
            return
        T0 = self.DeadStates['T0']
        P0 = self.DeadStates['P0']
        try:
            h0 = CP.PropsSI('H', 'T', T0, 'P', P0, self.fluid)
            s0 = CP.PropsSI('S', 'T', T0, 'P', P0, self.fluid)
        except Exception:
            object.__setattr__(self, 'ex', None)
            object.__setattr__(self, 'Ex', None)
            return

        ex = (self.H - h0) - T0 * (self.S - s0)
        Ex = self.Mass_flowrate * ex if self.Mass_flowrate else None
        object.__setattr__(self, 'ex', ex)
        object.__setattr__(self, 'Ex', Ex)

    # ------------------------------------------------------------------ #
    #  Therminol-66 polynomial calculation
    # ------------------------------------------------------------------ #
    @staticmethod
    def _t66_entropy(cp_coeffs, entropy_spec, T):
        """
        Specific entropy of an incompressible liquid, s(T) = ∫ Cp/T dT.

        For a quadratic Cp = a·T² + b·T + c this integrates to

            s(T) - s_ref = a/2·(T² - T_ref²) + b·(T - T_ref) + c·ln(T/T_ref)

        which is NOT a polynomial — the logarithmic term is why the original
        ``T66.json`` could not carry an ``entropy`` polynomial and why the
        key was missing altogether (raising ``KeyError: 'entropy'`` for every
        Therminol-66 state point up to v3.2.1).

        The coefficients are read from the ``specific_heat`` block, so the
        enthalpy, Cp and entropy correlations can never drift out of step.
        """
        T_ref = entropy_spec.get('T_ref', 273.15)
        s_ref = entropy_spec.get('s_ref', 0.0)

        # np.polyval order is highest power first.
        c = list(cp_coeffs)
        while len(c) < 3:
            c.insert(0, 0.0)
        a, b, c0 = c[-3], c[-2], c[-1]

        return (s_ref
                + 0.5 * a * (T ** 2 - T_ref ** 2)
                + b * (T - T_ref)
                + c0 * np.log(T / T_ref))

    def _calc_therminol66(self):
        t66_path = None
        for p in _T66_PATHS:
            if os.path.isfile(p):
                t66_path = p
                break
        if t66_path is None:
            raise FileNotFoundError(
                "T66.json not found. Looked in: " + ", ".join(_T66_PATHS)
            )

        with open(t66_path, 'r') as f:
            prop = json.load(f)

        data   = prop['fluids']['T66']
        h_mean = prop['h_normalization']['h_mean']
        h_std  = prop['h_normalization']['h_std']

        required = ['density', 'specific_heat', 'enthalpy', 'inverse_enthalpy']
        missing  = [k for k in required if k not in data]
        if missing:
            raise KeyError(
                f"T66.json is missing correlation block(s): {missing}. "
                f"Present: {sorted(data)}"
            )

        coeffs = {k: data[k]['coeffs'] for k in required}
        entropy_spec = data.get('entropy', {'type': 'cp_integral'})

        # Guard: disable hook during internal writes
        object.__setattr__(self, '_initialising', True)
        try:
            # --- resolve temperature ------------------------------------
            # Whichever of (T, H) the user supplied is authoritative and is
            # NOT overwritten: the enthalpy correlation and its inverse are
            # independent fits, so round-tripping would silently perturb the
            # user's input and break downstream energy balances.
            T_given = self.T is not None
            H_given = self.H is not None

            if T_given:
                pass
            elif H_given:
                H_norm = (self.H - h_mean) / h_std
                object.__setattr__(
                    self, 'T',
                    float(np.polyval(coeffs['inverse_enthalpy'], H_norm)))
            else:
                raise ValueError(
                    f"State '{self.StatePointName}': Therminol66 requires "
                    f"T or H as one of the two inputs."
                )

            T = self.T
            T_min = data.get('Tmin', -3.0) + 273.15
            T_max = data.get('Tmax', 350.0) + 273.15
            if not (T_min <= T <= T_max):
                warnings.warn(
                    f"State '{self.StatePointName}': T = {T - 273.15:.1f} °C "
                    f"is outside the fitted Therminol-66 range "
                    f"({T_min - 273.15:.0f} – {T_max - 273.15:.0f} °C). "
                    f"Properties are extrapolated and may be unreliable."
                )

            object.__setattr__(self, 'D',  float(np.polyval(coeffs['density'],      T)))
            object.__setattr__(self, 'Cp', float(np.polyval(coeffs['specific_heat'], T)))
            if not H_given:
                object.__setattr__(
                    self, 'H', float(np.polyval(coeffs['enthalpy'], T)))
            object.__setattr__(
                self, 'S',
                float(self._t66_entropy(coeffs['specific_heat'], entropy_spec, T)))

            # Therminol-66 is modelled as a single-phase incompressible
            # liquid: vapour quality is undefined.
            object.__setattr__(self, 'Q', None)
            object.__setattr__(self, 'phase', 'Incompressible liquid')

            # exergy
            if self.DeadStates:
                T0 = self.DeadStates['T0']
                h0 = float(np.polyval(coeffs['enthalpy'], T0))
                s0 = float(self._t66_entropy(
                    coeffs['specific_heat'], entropy_spec, T0))
                ex = (self.H - h0) - T0 * (self.S - s0)
                Ex = (self.Mass_flowrate * ex
                      if self.Mass_flowrate is not None else None)
                object.__setattr__(self, 'ex', ex)
                object.__setattr__(self, 'Ex', Ex)

            object.__setattr__(self, 'Solution_Status', True)

        finally:
            object.__setattr__(self, '_initialising', False)

    # ------------------------------------------------------------------ #
    #  public helper (used by Separator, etc.)
    # ------------------------------------------------------------------ #
    @staticmethod
    def _is_therminol66(fluid):
        """True when `fluid` names the built-in Therminol-66 correlation."""
        return isinstance(fluid, str) and fluid.strip().lower() in _T66_ALIASES

    def calculate_missing_properties(self):
        if self._is_therminol66(self.fluid):
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
                f"  Phase     : {self.phase}\n"
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
