> **Status:** Implemented on the local `post-review` branch (not pushed;
> `main` mirrors the version under review at SoftwareX and is untouched).
> Phases 1, 2, 4 and 5 landed as planned. Phase 3 ended up fully
> implemented rather than deferred: building the worked example surfaced
> that `HeatExchanger`'s `SimpleHEX` mode needed zero changes (it never
> touches the pinch machinery), and the general pinch-based solver
> (`double_pipe`/`Condenser`/`Evaporator`) was made mixture-aware too,
> including its 2- and 3-unknown cases. Two follow-up items outside this
> plan's original scope were also completed in the same branch: a
> pre-existing (non-mixture-specific) solver bug found while testing
> Phase 3, and a `work=` argument for `Turbine`/`Compressor`/`Pump`
> requested separately. Phase 6 (Kalina cycle) remains future work.

# Integrate NH3-H2O absorption cycle support into ThermoSim

## Context

`NH3H2O/nh3h2o.py` is a standalone, validated (22/22 checks, energy balance
closing to 6.6e-12% in `NH3H2O/demo_cycle.py`) NH3-H2O Gibbs-excess property
library the user wrote separately. They want to fold it into ThermoSim so
absorption-cycle and (later) Kalina-cycle models can be built with the same
`Model` / `Prop` / component API used for every other cycle in the package,
instead of hand-writing balance equations as `demo_cycle.py` currently does.

Because `main` mirrors the version currently under review at SoftwareX and
must not change, all of this work happens on the local-only `post-review`
branch (already created, currently 1 commit ahead of `main`). Nothing here
gets pushed.

Confirmed direction: build the **single-effect absorption chiller first**
(directly validatable against `demo_cycle.py`'s numbers); keep composition
**local to each state point** rather than extending the generic mass-flow
propagation graph; and — per the user's explicit steer — model the cycle by
**reusing existing generic components** wherever they already fit, instead
of writing bespoke `Absorber`/`Desorber` classes:

- **Absorber** = `Mixer` (weak solution + refrigerant vapor → mixed stream)
  followed by a `HeatExchanger` acting as a condenser (cools/subcools the
  mixed stream down to the strong-solution state, rejecting `Q_abs`).
- **Generator** = a `HeatExchanger` (heats the strong solution into the
  two-phase dome) followed by `Separator` used as-is.
- **Evaporator** / **Condenser** (refrigerant side) = `HeatExchanger`,
  used as-is once it's mixture-aware.

Kalina-cycle turbine/recuperator work stays an explicit follow-up phase.

## Architecture: composition rides on the fluid identifier

`ThermoSim/state.py`'s `Prop` class already special-cases one non-CoolProp
fluid — Therminol-66 — by string-matching `fluid` against `_T66_ALIASES`
(`state.py:44`) and branching to `_calc_therminol66()` instead of
`_calc_coolprop()`. NH3-H2O gets the same dispatch pattern, but also needs a
**composition**, which single-fluid CoolProp lookups never carried.

Rather than adding a `z` parameter every call site must remember to thread
through, `z` is bundled into the fluid identifier itself: a new small class
`NH3H2OFluid(z)` (hashable, `__repr__`-friendly) is what gets passed as
`fluid`. This matters because every existing component — `Pump`, `Turbine`,
`Compressor`, `Expansion_valve`, `Mixer` (outlet fluid is whatever the user
pre-declared on that point, `simple_components.py:277`), and the
`heat_exchangers.py` helpers `_h_from_T`/`_h_from_H`/`_T_sat` — already does
`self.Model.Prop(<some existing fluid>, P=..., S=..., ...)`, forwarding an
object it never inspects. With composition baked into that object, these
components work with the mixture with **no changes to their own logic**.

Two components *do* need small, targeted changes because they currently
assume "same fluid in, same fluid out" in a way that's only true for a
single species — both detailed below.

## Phase 1 — Mixture state layer

- Vendor `NH3H2O/nh3h2o.py` into the package as `ThermoSim/nh3h2o.py`
  (self-contained — only `numpy`/`scipy`, no CoolProp dependency). Keep
  `NH3H2O/README_nh3h2o.md`'s validity notes (230-500 K, 0.2-50 bar;
  non-uniqueness warnings for `(T,H)`/`(P,V)`) as the authoritative
  reference, folded into `docs/`.
- Add `NH3H2OFluid` (new small class, e.g. in `ThermoSim/nh3h2o.py`)
  wrapping a composition `z` (ammonia mass fraction), with `__repr__` like
  `NH3H2O(z=0.4000)`, `__eq__`/`__hash__` by `(name, z)`.
- In `ThermoSim/state.py`:
  - Add `Prop._is_nh3h2o(fluid)` (same shape as `_is_therminol66`,
    `state.py:483`).
  - Add `_calc_nh3h2o()`, following the `_calc_therminol66()` pattern
    (`state.py:385`): convert the two known SI inputs to `nh3h2o.py`'s units
    (Pa→bar for P; J/kg→kJ/kg for H; J/kg/K→kJ/kg/K for S), call
    `nh3h2o.state(z=fluid.z, **converted)`, convert back to SI, and populate
    `P, T, H, S, Q, D` plus new attributes `z, x, y` (equilibrium liquid/
    vapor compositions — `None` outside the two-phase dome, and equal to
    `z` trivially in the single-phase region).
  - Wire `_calc_exergy()` (`state.py:330`) to also handle the mixture case
    (dead-state `h0, s0` via `nh3h2o.state(z=fluid.z, T=T0, P=P0)`).
  - Map `nh3h2o`'s `AmmoniaWaterError` (the genuinely non-unique `(T,H)`/
    `(P,V)` inputs) into the same `ValueError` wrapping style
    `_calc_coolprop()` already uses.
  - `calculate_missing_properties()` (`state.py:488`) gains this third
    branch alongside CoolProp / Therminol-66.
- Export `NH3H2OFluid` from `ThermoSim/__init__.py` so a user writes
  `M.add_point(NH3H2OFluid(0.4), 'S1', P=..., T=...)`.

## Phase 2 — The one required fix: `Separator` must read equilibrium x/y

`Separator.Cal()` (`simple_components.py:531-612`) builds both outlets by
reusing the inlet's fluid unchanged: `self.Model.Prop(self.In.fluid, P=P,
Q=1.0, ...)` for vapor, `Q=0.0` for liquid. For a **pure** fluid this is
correct — saturated vapor/liquid at pressure P is one well-defined state
regardless of the actual inlet quality. For the **mixture** it is not:
`nh3h2o.state(z=z_feed, P=P, Q=1.0)` is the *dew point of the feed
composition* at P — a different point on the envelope than the vapor
phase actually coexisting with the liquid at the real flash temperature.
The physically correct outlet compositions are the equilibrium `x`/`y` the
flash already computed on the inlet `Prop` (Phase 1 adds these).

So `Separator` needs one small, additive change: when `self.In.fluid` is an
`NH3H2OFluid`, build the two outlets with `NH3H2OFluid(self.In.x)` /
`NH3H2OFluid(self.In.y)` instead of `self.In.fluid` verbatim, and flash
them at `(P, T_actual)` rather than the `Q=0`/`Q=1` boundary shortcut
(equivalently: `(P, self.In.T)`, since liquid and vapor are in equilibrium
at the same T). Everything else in `Separator` — the tiered P/mass-flow
logic, the energy-balance sanity check — is untouched and reused as-is,
matching the spirit of the user's direction.

`Mixer` needs no code change to work correctly for the absorber-mixing
step — its outlet composition is always whatever the user pre-declared on
the `Out_state` point (`self.Out.fluid`, set at `add_point()` time), and
`demo_cycle.py` already computes that target composition externally via a
species-balance solve (`z_s = brentq(...)`, `demo_cycle.py:20`) before any
state is built — the worked example in Phase 4 does the same. One
*optional* hardening worth adding: a warning if `Σ(m_i·z_i)/Σm_i` (computed
from the inlets) deviates from the user-declared outlet `z`, mirroring the
enthalpy-consistency check `Mixer` already does (`simple_components.py:
268-274`) — catches a mismatched target composition instead of silently
building a state that violates species conservation.

## Phase 3 — Mixture support in `HeatExchanger`

This is the piece the "Absorber = Mixer + condenser-HEX" /
"Generator = HEX + Separator" / "Evaporator = HEX" design depends on, and
the biggest remaining engineering item: `heat_exchangers.py`'s pinch-point
solver has two helpers that call `CP.PropsSI` **directly**, bypassing
`Prop` entirely — `_safe_T_min` (`heat_exchangers.py:48`) and
`_phase_boundary_alphas` (`:175`) — plus (to be confirmed while reading the
~1800-line `HeatExchanger.Cal()` in full during implementation) likely
further direct CoolProp calls in the pinch-profile temperature scan itself.
These assume a **flat saturation temperature** at a given P, which is true
for a pure fluid but not for a binary mixture (NH3-H2O has a temperature
*glide* across its two-phase dome at fixed P).

Proposed staged approach, to avoid blocking the first validated cycle on
the hardest part of the codebase:

1. First pass: use `HeatExchanger` only in configurations where the
   mixture-side inlet **and** outlet states are already fully determined by
   the flash relations (exactly how every mixture state in `demo_cycle.py`
   is obtained — e.g. the strong solution enters the generator at the pump
   outlet state and leaves at the known two-phase flash state at `T_gen`).
   In that configuration the HEX only needs to solve the utility-side
   flow/duty from `Q = m·Δh` on the known side — no pinch-profile search
   over the mixture's glide is required yet, so `_safe_T_min` and
   `_phase_boundary_alphas` are not on the critical path for a first
   working, validated cycle.
2. Extend `_safe_T_min`, `_phase_boundary_alphas`, and the internal
   pinch-scan to dispatch to `nh3h2o` for `NH3H2OFluid` streams (computing
   the actual glide curve instead of a single saturation point), so
   `HeatExchanger` can later auto-solve harder configurations (unknown
   mixture-side terminal temperature, multi-unknown cases) the same way it
   already does for every CoolProp fluid. This is real, separate work and
   is called out explicitly rather than folded silently into step 1.

## Phase 4 — Worked example, validated against `demo_cycle.py`

Add `examples/nh3h2o_absorption_chiller.py` (naming convention matching
`validation/fig_chena.py`) that builds the single-effect chiller from
`ThermodynamicModel` + `Prop`/`NH3H2OFluid` + `Pump` + `Expansion_valve` +
`Mixer` + `HeatExchanger` + `Separator`, following `demo_cycle.py`'s state
sequence exactly, and asserts `COP`, `Q_gen`, `Q_cond`, `Q_evap`, `Q_abs`
match `demo_cycle.py`'s printed reference values (already captured in
`NH3H2O/validation_output.txt`) to a tight tolerance.

The rectifier stays a fixed-purity assumption (`y_r = 0.995`, as
`demo_cycle.py` itself declares rather than rigorously solving a
multi-stage column) — out of scope for this pass, noted as future work
alongside Kalina.

## Phase 5 — Tests

- `tests/test_nh3h2o_state.py`: round-trip checks mirroring
  `NH3H2O/validate_lib.py`'s existing 22 checks, run through `Prop`/
  `NH3H2OFluid` rather than the bare library, confirming the SI
  unit-conversion layer is exact.
- `tests/test_absorption_cycle.py`: builds the Phase 4 example through
  `Model.Solve()` and asserts mass/species/energy balances close to the
  same tolerance `demo_cycle.py` reports.
- Run the existing suite (`pytest tests/`) to confirm the new `Prop`
  dispatch branch and the `Separator` change don't regress any existing
  single-fluid cycle — the key safety check for touching shared code paths
  (`Prop.__setattr__`'s auto-recalc hook, `Separator.Cal()`).

## Phase 6 (future, out of scope) — Kalina cycle

Same `NH3H2OFluid`/state-layer foundation; adds turbine expansion of
superheated mixture vapor (`Turbine`, unmodified — same fluid-forwarding
argument as `Pump`) and recuperators (another `HeatExchanger` mixture
use-case, benefiting from Phase 3 step 2). Flag up front: Kalina cycles
often run above the Ziegler-Trepp validity envelope `nh3h2o.py` documents
(230-500 K, ≤50 bar, accuracy degrading sharply above that) — needs a
validity check/warning at minimum, possibly a different high-pressure
correlation later (Ibrahim & Klein, or Tillner-Roth & Friend, both cited in
`NH3H2O/README_nh3h2o.md`).

## Verification

- Run `python examples/nh3h2o_absorption_chiller.py`, diff its printed
  COP/duties against `NH3H2O/demo_cycle.py` / `NH3H2O/validation_output.txt`.
- Run `tests/test_nh3h2o_state.py` and `tests/test_absorption_cycle.py`.
- Run `pytest tests/` in full for regression safety.
- Commit each phase separately on `post-review`, never touching `main`.
