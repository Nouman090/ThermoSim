# Contributing to ThermoSim

Thanks for your interest. Bug reports, questions and pull requests are all
welcome.

## Reporting a bug

Open an issue at
<https://github.com/Nouman090/ThermoSim/issues> and include:

1. **A minimal script that reproduces it.** The smallest model that still
   shows the problem — usually a handful of state points and one or two
   components.
2. **What you expected and what happened**, including the full traceback if
   there was one.
3. **Your environment**:

   ```bash
   python -c "import sys, ThermoSim, CoolProp; \
   print(sys.version.split()[0], ThermoSim.__version__, CoolProp.__version__)"
   ```

Physically wrong answers matter as much as crashes. If a cycle solves but
the numbers look implausible, that is worth reporting — several defects
found during validation produced no error at all.

## Asking a question

Open an issue and label it a question. Questions about how to model a
particular cycle are welcome; they often reveal gaps in the documentation.

## Contributing code

1. Fork the repository and create a branch off `main`.
2. Make the change.
3. Run the test suite:

   ```bash
   python -m pytest tests/ -q          # expect 57 passed
   ```

4. If the change touches the solver, also run the validation scripts:

   ```bash
   cd validation
   python validate_moran.py
   python validate_moran_86.py
   python validate_moran_gas.py
   python validate_moran_refrig.py
   python validate_chena.py
   ```

   These compare against published worked examples and measured plant data.
   Any change in the reported deviations should be understood and explained
   before the pull request is opened.

5. Add a test for the behaviour you changed.
6. Open a pull request describing what changed and why.

Continuous integration runs the test suite on Python 3.9 through 3.13 on
Linux for every pull request.

## Development install

```bash
git clone https://github.com/Nouman090/ThermoSim.git
cd ThermoSim
pip install -e .
pip install pytest
python -m pytest tests/ -q
```

## Style

- Follow the conventions already in the file you are editing.
- All internal quantities are SI: Pa, K, J/kg, J/(kg·K), kg/s, W, J. Unit
  conversion belongs in printing and plotting, not in the solver.
- Prefer a clear error message over a silent fallback. If a specification
  cannot be solved, say which combination was supplied and why it does not
  close.
- Docstrings explain *why*, not just *what*. The reasoning behind a
  non-obvious choice is the part a future reader cannot reconstruct.

## Scope

ThermoSim models thermodynamic cycles at the state-point level. Some things
are deliberately outside its scope:

- **Combustion chemistry.** Reactions are represented as heat input.
- **Heat exchanger geometry.** `UA` is reported; sizing is not.
- **Transient behaviour.** All analysis is steady-state.

Proposals that extend these boundaries are worth discussing in an issue
first.

## Licence

Contributions are accepted under the MIT Licence, the same terms as the
rest of the project.
