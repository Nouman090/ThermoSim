
# Changelog

## [3.0.0] - 2026-03-25

### Added
- Modular package structure (`thermocycle/` package with separate files)
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
