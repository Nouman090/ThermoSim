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
import numpy as np
import CoolProp.CoolProp as CP

from . import config                       # reads dead_states from config.py
from . import nh3h2o as _aw
from .nh3h2o import NH3H2OFluid, AmmoniaWaterError as _AmmoniaWaterError

# Attributes whose assignment triggers a full thermodynamic re-solve.
_THERMO_INPUTS = {'P', 'T', 'H', 'S', 'Q', 'D'}

# Short names accepted for Therminol 66.  Up to v3.5.0 ThermoSim carried its
# own polynomial fit (T66.json) for these; CoolProp's incompressible backend
# provides Therminol 66 directly, so every alias now resolves to it.
# Matching is case-insensitive.
_FLUID_ALIASES = {
    'therminol66':  'INCOMP::T66',
    'therminol-66': 'INCOMP::T66',
    't66':          'INCOMP::T66',
    'incomp::t66':  'INCOMP::T66',
}


def canonical_fluid(fluid):
    """Return the CoolProp name for `fluid`, resolving ThermoSim aliases.

    Non-string fluids (an NH3H2OFluid) and names without an alias are
    returned unchanged.
    """
    if isinstance(fluid, str):
        return _FLUID_ALIASES.get(fluid.strip().lower(), fluid)
    return fluid


def is_incompressible(fluid):
    """True for CoolProp incompressible-backend fluids (``INCOMP::...``).

    These are single-phase liquids: CoolProp has no saturation states or
    vapour quality for them.
    """
    return (isinstance(fluid, str)
            and canonical_fluid(fluid).strip().upper().startswith('INCOMP::'))

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

        self.fluid = canonical_fluid(fluid)
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
        # NH3-H2O mixture composition: overall (z) and, inside the
        # two-phase dome, the coexisting liquid (x) / vapor (y) compositions.
        # None for every other fluid type.
        self.z = fluid.z if isinstance(fluid, NH3H2OFluid) else None
        self.x = self.y = None

        # set the properties we already know
        for key, value in self.properties.items():
            setattr(self, key, value)

        # --- calculate the rest -----------------------------------------
        if len(self.properties) == 2:
            if self._is_nh3h2o(self.fluid):
                self._calc_nh3h2o()
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
                if self._is_nh3h2o(self.__dict__.get('fluid')):
                    self._calc_exergy_nh3h2o()
                else:
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
        incomp = is_incompressible(self.fluid)

        # --- handle near-saturation for P-T inputs ---------------------
        # CoolProp refuses a (T, P) flash within ~1e-4 % of the saturation
        # pressure, so inputs within a small relative band of T_sat are
        # snapped onto the saturation line.  The snap is SIDE-AWARE: a T
        # above T_sat is saturated vapour (Q=1), a T below it (or equal to
        # it within round-off) is saturated liquid (Q=0).  Up to v3.5.0 it
        # always chose Q=0, turning a slightly superheated suction state
        # into liquid.  For a zeotropic mixture, T_sat (Q=0) is the bubble
        # point and a T just above it is two-phase with tiny quality, so
        # Q=1 is chosen only when T is nearer the dew point than the bubble
        # point (always true for a pure fluid, where they coincide).
        if {prop1, prop2} == {'P', 'T'} and not incomp:
            P = self.properties['P']
            T = self.properties['T']
            try:
                T_sat = CP.PropsSI('T', 'P', P, 'Q', 0, self.fluid)
                if abs(T - T_sat) / T_sat < 1e-4:
                    Q_snap = 0.0
                    # Round-off band: T == T_sat(P) keeps the historical Q=0.
                    if (T - T_sat) / T_sat > 1e-7:
                        T_dew = CP.PropsSI('T', 'P', P, 'Q', 1, self.fluid)
                        if T - T_sat > T_dew - T:
                            Q_snap = 1.0
                    side = ("saturated vapour (Q=1)" if Q_snap == 1.0
                            else "saturated liquid (Q=0)")
                    warnings.warn(
                        f"State '{self.StatePointName}' is near saturation "
                        f"(T - T_sat = {T - T_sat:+.4f} K). Assuming {side}."
                    )
                    prop1, val1 = 'P', P
                    prop2, val2 = 'Q', Q_snap
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
            # Incompressible liquids have no vapour quality (CoolProp raises).
            if self.Q  is None and not incomp:
                object.__setattr__(self, 'Q',  CP.PropsSI('Q', prop1, val1, prop2, val2, self.fluid))
            if self.D  is None: object.__setattr__(self, 'D',  CP.PropsSI('D', prop1, val1, prop2, val2, self.fluid))
            if self.Cp is None: object.__setattr__(self, 'Cp', CP.PropsSI('C', prop1, val1, prop2, val2, self.fluid))

            if incomp:
                object.__setattr__(self, 'Q', None)
                object.__setattr__(self, 'phase', 'Incompressible liquid')
            else:
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
    #  public helper (used by Separator, etc.)
    # ------------------------------------------------------------------ #
    @staticmethod
    def _is_nh3h2o(fluid):
        """True when `fluid` is an NH3H2OFluid (binary mixture) tag."""
        return isinstance(fluid, NH3H2OFluid)

    def calculate_missing_properties(self):
        if self._is_nh3h2o(self.fluid):
            self._calc_nh3h2o()
        else:
            self._calc_coolprop()

    # ------------------------------------------------------------------ #
    #  NH3-H2O mixture calculation
    # ------------------------------------------------------------------ #
    def _calc_nh3h2o(self):
        """
        Flash an NH3-H2O state via ``nh3h2o.py`` and translate between
        ThermoSim's SI unit convention and the library's own (P in bar,
        H/S/U in kJ/kg[/K]). T, Q and density (v = 1/D, both already SI
        m3/kg <-> kg/m3) need no conversion.

        Composition is carried on ``self.fluid`` (an ``NH3H2OFluid``), not
        on ``self.properties`` -- the two thermodynamic inputs here are
        exactly the two keys already validated by ``_VALID_PROPS``.
        """
        prop1, prop2 = list(self.properties.keys())
        val1, val2 = list(self.properties.values())

        _TO_AW = {
            'P': lambda v: v / 1e5,
            'H': lambda v: v / 1e3,
            'S': lambda v: v / 1e3,
            'T': lambda v: v,
            'Q': lambda v: v,
            'D': lambda v: v,
        }
        kwargs = {prop1: _TO_AW[prop1](val1), prop2: _TO_AW[prop2](val2)}

        object.__setattr__(self, '_initialising', True)
        try:
            try:
                st = _aw.state(z=self.fluid.z, **kwargs)
            except _AmmoniaWaterError as e:
                raise ValueError(
                    f"NH3-H2O flash error for state '{self.StatePointName}': {e}"
                )

            object.__setattr__(self, 'P', st.P * 1e5)
            object.__setattr__(self, 'T', st.T)
            object.__setattr__(self, 'H', st.h * 1e3)
            object.__setattr__(self, 'S', st.s * 1e3)
            object.__setattr__(self, 'D', 1.0 / st.v)
            object.__setattr__(self, 'z', st.z)

            if st.phase == 'two-phase':
                object.__setattr__(self, 'Q', st.Q)
                object.__setattr__(self, 'x', st.x)
                object.__setattr__(self, 'y', st.y)
                object.__setattr__(
                    self, 'phase', f"Two-phase (Q={st.Q:.4f})")
            else:
                object.__setattr__(self, 'Q', None)
                object.__setattr__(self, 'x', st.z)
                object.__setattr__(self, 'y', st.z)
                label = 'Liquid solution' if st.phase == 'liquid' else 'Vapor'
                object.__setattr__(self, 'phase', label)

            # Cp is not provided by nh3h2o.py (no cp field on State); left
            # None, matching every other attribute this library doesn't
            # compute rather than silently guessing at it.
            object.__setattr__(self, 'Cp', None)

            self._calc_exergy_nh3h2o()

            object.__setattr__(self, 'Solution_Status', True)
        finally:
            object.__setattr__(self, '_initialising', False)

    def _calc_exergy_nh3h2o(self):
        """Specific/total flow exergy for an NH3-H2O state, mirroring
        ``_calc_exergy`` but flashing the dead state through ``nh3h2o.py``
        (with the SAME overall composition) instead of CoolProp."""
        if not self.DeadStates or self.H is None or self.S is None:
            object.__setattr__(self, 'ex', None)
            object.__setattr__(self, 'Ex', None)
            return
        T0 = self.DeadStates['T0']
        P0 = self.DeadStates['P0']
        try:
            st0 = _aw.state(z=self.fluid.z, T=T0, P=P0 / 1e5, warn=False)
            h0, s0 = st0.h * 1e3, st0.s * 1e3
        except _AmmoniaWaterError:
            object.__setattr__(self, 'ex', None)
            object.__setattr__(self, 'Ex', None)
            return

        ex = (self.H - h0) - T0 * (self.S - s0)
        Ex = self.Mass_flowrate * ex if self.Mass_flowrate else None
        object.__setattr__(self, 'ex', ex)
        object.__setattr__(self, 'Ex', Ex)

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
