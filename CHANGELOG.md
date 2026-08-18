
# Changelog

## [3.2.2] - 2026-08-07

### Fixed — packaging
- `tabulate` added to the declared dependencies. It is imported at module
  scope in `model.py`, so a clean `pip install ThermoSim` followed by
  `import ThermoSim` previously raised `ModuleNotFoundError`.
- Placeholder author name/email replaced in `pyproject.toml`.
- `heat_exchangers_old.py` (dead, never imported) removed from the package.
- `.gitignore` added; it excludes private keys, `.env` files and build output.
- `README.md` no longer contains unresolved Git merge-conflict markers.

### Fixed — state points (`Prop`)
- `Q` is now ALWAYS numeric (0-1) or `None`. The human-readable phase label
  moved to the new `Prop.phase` attribute. Previously `_classify_phase()`
  overwrote `Q` with strings like `"Superheated (46.97 K)"`, which crashed
  `Separator` and made any arithmetic on `Q` unsafe.
- Assigning the FIRST property to an under-defined point no longer raises
  `ValueError: not enough values to unpack`; a solve is attempted only once
  two independent inputs are present.
- Therminol-66 works again. `T66.json` was missing the `entropy` block that
  `_calc_therminol66()` required, so every Therminol-66 state raised
  `KeyError: 'entropy'`. Entropy is now integrated analytically from the Cp
  polynomial (s = ∫Cp/T dT, which contains a log term and therefore cannot
  be a polynomial). Whichever of (T, H) the user supplies is preserved
  rather than round-tripped through the independent inverse fit, the fitted
  temperature range is range-checked, and `'T66'` is accepted as an alias.

### Fixed — turbomachinery
- `Turbine`, `Compressor` and `Pump` write the resolved mass flow back to
  their ports again. Since 3.1.1 the flow was kept local and left to the
  propagation graph, which is only built by `Solve()` — so the documented
  `Calculate=True` workflow left every downstream point at
  `Mass_flowrate=None` and the next heat exchanger could not solve.
- Mass-flow resolution uses an explicit `is not None` test, so a legitimate
  0 kg/s is no longer treated as missing.

### Fixed — heat exchangers
- Mass flows are resolved BEFORE the SimpleHEX branch. Previously that branch
  returned early, so `Hot/Cold_Mass_flowrate` were still `None` and the
  solver took the wrong path.
- `_solve_simple_hex()` rewritten: each side is now closed independently
  against a shared duty, so a two-sided unit no longer solves only one side;
  a fully-determined unit is idempotent (calling `Cal()` twice used to raise
  "insufficient inputs", which made every `Solve()` iteration report phantom
  failures); and an over-specified `Q` is checked rather than ignored.
- The NTU-effectiveness solver is reachable. It was the last `elif` in a
  chain that had already consumed n = 0, 1, 2 and 3, so the `effectiveness`
  argument silently did nothing.
- A failed 3-unknown solve no longer corrupts the model: the anchor is
  validated against the solvable-case table first and rolled back on failure,
  and both valid anchors are tried before giving up.
- Anchor rules no longer compare `None` temperatures (`TypeError`), and the
  `H_co → H_hi` guard tests the correct temperature pair; it previously
  reused the `H_ho → H_ci` condition verbatim.
- The `{H_co, mc}` retry passes the narrowed bracket instead of re-using the
  bracket that had just failed, and cold-side enthalpies are no longer
  bounded by hot-side enthalpies (meaningless across different fluids).
- `exclude_idx` masking now covers the final profile index, which is exactly
  the anchored point it was meant to exclude.
- `_compute_outputs()` uses `is not None` instead of truthiness, so a valid
  H = 0 or Q = 0 no longer silently skips the duty and UA calculations.
- Exergy destruction is computed per connected leg, so one-sided units report
  a value instead of falling back to "Not Calculated" and printing to stdout.
- Solvable-case counts in the module docstring corrected to 10 (2-unknown)
  and 8 (3-unknown), matching the logic table and the CHANGELOG. The list is
  now a single constant, `HeatExchanger._SOLVABLE_2_UNKNOWN`.
- `HeatExchanger` and `TES` call `super().__init__()`; `TES` reuses the
  shared mass-flow resolver instead of comparing floats with `==`.

### Fixed — model
- `ModelSummary()` skips unsolved components with a warning instead of
  raising `TypeError` on `Q_in += None`.
- `Efficiency` is `nan` rather than the string `"Not Applicable"` for
  refrigeration cycles, and `COP_R` / `COP_HP` are exposed as attributes.
- `Solve()` convergence uses a relative tolerance (new `tol` argument). The
  old absolute 1e-6 J/kg was ~3e-13 relative on typical enthalpies.
- `Solve()` and the flow network are quiet unless `verbose=True`.
- `custom_order` defaults to `None` instead of a mutable list.
- `Point_print()` no longer advertises private attributes.

### Fixed — sensitivity analysis
- `compute_sensitivity_indices()` returned a meaningless `sensitivity_index`.
  It normalised BOTH axes to [0, 1] before taking the gradient, so any
  monotonic near-linear sweep gave ~1.0000 regardless of the parameter — the
  metric could not rank parameters against each other, which is the one thing
  the wiki instructs users to do with it. It is now the mean **elasticity**,
  S = (dY/dX)·(X/Y): dimensionless, so parameters in kelvin and pascals stay
  comparable. Validation on a simple Rankine cycle now gives eta_t → 1.01
  (cycle efficiency is nearly proportional to turbine efficiency, as
  expected), T1 → 0.27, P_cond → 0.07. The raw slope is also returned, as
  the new `mean_gradient` key.
- `ModelSummary()` accepts `verbose=True` (default). `SensitivityAnalyzer`
  calls it with `verbose=False`, so a sweep no longer prints a full
  energy + exergy table for every single evaluation point.

### Fixed — components, plotting, analysis
- `Splitter` treats a known inlet flow as authoritative instead of averaging
  it with stale outlet flows. Re-running a splitter with refined fractions
  (as the regenerative Rankine example does) previously broke mass
  conservation downstream.
- Every component `__str__` returns "not yet solved" instead of raising.
- `plot_energy_summary()` no longer raises `NameError` on `Compressor` in its
  import-failure fallback, and tolerates unsolved components.
- `'Water'` removed from `INCOMPRESSIBLE_FLUIDS`; it is an ordinary
  compressible CoolProp fluid and was silently losing its saturation dome.
- `plot_hex_profile()` reports a clear error instead of dividing by `None`.
- `compute_sensitivity_indices()` guards against a zero-width parameter range.
- Bare `except:` replaced with `except Exception:` throughout (22 sites).
- Corrected `thermocycle.*` import paths in docstrings to `ThermoSim.*`.

## [3.2.0] - 2026-04-20
- Solver is now more robust.
- User can give the order of solving components.
- Now HEX can solve 10 cases with 2 unknows and 8 cases with 3 unknwons
- Mass flowrate now can automatically propgate through the branches. 

## [3.1.1] - 2026-03-31
 - Fix the bug in testing.
## [3.1.0] - 2026-03-30
- Energy-balance consistency check when all four states are known.
- Mass-flow written back after PPT/brentq solvers (was silently discarded).
- Upfront problem classifier (_classify_problem) catches zero-info,
        both-m̊-unknown, and over-specified inputs before any computation.
- Silent PPT mutation (self.PPT -= 0.001) replaced with a proper
        ValueError so pinch violations are never hidden.
- Phase-change guard in _solve_effectiveness(); falls back to
        enthalpy-based Q when either stream is two-phase (Cp is ill-defined).
- TES over-specified branch added: CI.H, CO.H, and Capacity all known
        → consistency check instead of silent fall-through.
- warnings.warn(RuntimeWarning) replaces bare print() in solver errors;
        Solution_Status set False and Cal() returns cleanly on non-fatal failures.
- SimpleHEX two-sided consistency path: when both cold and hot sides are
        fully given, Q from each side is compared and averaged if consistent.
-  _resolve_side() raises immediately when both mass-flows are None
        (was a silent None that propagated to arithmetic errors later).
- exergy sign and base-class helper preserved from v1.x bugfix baseline.
- All three Cal() methods now handle a "both-known" branch that validates
    isentropic consistency and emits a warning when the supplied enthalpy
    pair violates the declared efficiency.
- Pump.Cal() Incompressible branch now handles reverse (In.H is None)
    and both-known sub-cases with proper guards.
- fsolve calls use physically-motivated initial guesses (rather than
    the outlet enthalpy) and check convergence with full_output=True.
- Guard added in the Incompressible path for In.D being None.
- Duplicate 'from scipy.optimize import fsolve' inside Compressor.Cal()
    removed.
- All three __str__ methods guard against Solution_Status=False (i.e.
    Cal() was never called) so printing an unsolved component does not
    raise AttributeError on self.In / self.Out.

- CSV Export Functionality

Added CSV export capability to all plotting methods

    - save_csv parameter (default: False) - enables CSV export
    - csv_path parameter (optional) - custom file path, auto-generates if not provided
    - Auto-generated filenames include timestamp (e.g., Ts_diagram_20240115_143022.csv)
New utility methods:

    - export_all_points() - Export all thermodynamic state points to CSV
    - export_all_components() - Export all component performance data to CSV
    - _generate_csv_filename() - Auto-generate timestamped filenames
    - _save_dataframe_to_csv() - Centralized CSV saving with error handling
    Enhanced Data Structure

CSV files now include:
    - Original SI units (J/kg, Pa, K)
    - Converted engineering units (kJ/kg, bar, °C)
    - Point/component identification
    - Calculated derivatives (percentages, ratios)

## [3.0.0] - 2026-03-25

### Added
- Modular package structure (`ThermoSIm/` package with separate files)
- `Compressor` component for gas cycles (e.g., Brayton)
- `CyclePlotter` class with T-s, P-h, h-s diagrams, saturation domes
- Exergy destruction bar chart and pie chart
- Heat exchanger temperature profile plots
- Energy flow summary visualisation
- `SensitivityAnalyzer` with single sweep, double sweep, and multi-output
- `save_model()` and `load_model()` for JSON serialisation
- Abstract `Component` base class with shared helpers
- 48+ unit tests (pytest)
- `.gitignore`, `README.md`, `CHANGELOG.md`, `setup.py`

### Fixed
- `__init__` was misspelled as `init` (constructor never called)
- Exergy sign error in HeatExchanger and TES: was `Out_hot - Out_cold`, now `Out_hot + Out_cold`
- `ModelSummary()` crashed on second call (DataFrame append issue)
- Type checking used fragile string comparison instead of `isinstance()`
- `Solve()` silently swallowed all errors
- Multiple typos in error messages and `__str__` methods

### Changed
- Dead state stored in `config.py` module instead of class variable
- Mass-flow resolution logic extracted to base class
- Visualisation completely separated from computation
- Iterative solver now checks convergence

## [2.0.0] - Previous version
- Single-file module

## [1.0.0] - Initial release
- Basic functionality
