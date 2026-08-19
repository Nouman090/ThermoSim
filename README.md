[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22016326.svg)](https://doi.org/10.5281/zenodo.22016326)
![tests](https://github.com/Nouman090/ThermoSim/actions/workflows/tests.yml/badge.svg)
# ThermoSim

![ThermoSim logo](https://github.com/Nouman090/ThermoSim/blob/main/docs/ThermoSim%20Logo%202.jpg?raw=true)

**ThermoSim** is a Python package for modelling and analysing thermodynamic
power cycles, built on [CoolProp](http://www.coolprop.org/) for fluid
properties.

---

## Features

- **State-point calculation** from any two independent properties
  (P, T, H, S, Q, D) for water, refrigerants, air and every other CoolProp
  fluid, plus a built-in Therminol-66 correlation
- **Components**: Turbine, Pump, Compressor, Heat Exchanger, TES, Mixer,
  Splitter, Separator, Expansion Valve, Pipe, Source, Sink
- **Pinch analysis** for heat exchangers, with a generalised solver that
  closes up to three unknowns from
  `{H_hi, H_ho, H_ci, H_co, mh, mc}`
- **Automatic mass-flow propagation** through branches, splitters and mixers
- **Exergy analysis** against a configurable dead state
- **Cycle plots**: T-s, P-h and h-s diagrams with saturation domes, heat
  exchanger temperature profiles, exergy bar/pie charts, energy summaries
- **Sensitivity analysis**: single-parameter, multi-output and 2D contour
  sweeps, with CSV export
- **Save / load** model state to JSON
- **57 unit tests** (pytest)

---

## Installation

```bash
pip install ThermoSim
```

To install from source:

```bash
git clone https://github.com/Nouman090/ThermoSim.git
cd ThermoSim
pip install -e .
```

Optional extras: `pip install "ThermoSim[all]"` adds `pymoo` (optimisation)
and `pytest` (testing).

---

## Quick example

```python
from ThermoSim import ThermodynamicModel, Turbine, Pump, HeatExchanger
from ThermoSim.plotting import CyclePlotter

Model = ThermodynamicModel()
Model.set_dead_state()

Model.add_point('water', '1', P=8e6,     T=753.15, Mass_flowrate=1)
Model.add_point('water', '2', P=0.008e6)
Model.add_point('water', '3', P=0.008e6, Q=0)
Model.add_point('water', '4', P=8e6)

Turbine(Model, 'Turbine', '1', '2', n_isen=0.85, Calculate=True)
HeatExchanger(Model, 'Condenser', PPT=5, HEX_type='SimpleHEX',
              HeatAdded=False, Hot_In_state='2', Hot_Out_state='3',
              Cold_In_state=None, Cold_Out_state=None, Calculate=True)
Pump(Model, 'Pump', '3', '4', n_isen=1.0, Calculate=True)
HeatExchanger(Model, 'Boiler', PPT=5, HEX_type='SimpleHEX',
              HeatAdded=True, Hot_In_state=None, Hot_Out_state=None,
              Cold_In_state='4', Cold_Out_state='1', Calculate=True)

print(Model)

plotter = CyclePlotter(Model)
plotter.plot_Ts_diagram(['1', '2', '3', '4', '1'])
```

For cycles with recycle loops or branches, build the components without
`Calculate=True` and let the iterative solver handle the ordering:

```python
Model.Solve(verbose=True)
Model.ModelSummary()
```

---

## Running the tests

```bash
pytest tests/ -v
```

---

## Who this is for

- **Engineering students** studying thermodynamics, energy systems and heat
  transfer, who want to simulate real energy systems rather than only work
  through closed-form textbook problems.
- **Researchers** modelling, optimising and analysing thermodynamic and
  energy-efficiency problems.
- **Energy system designers** working on power plants, heat exchanger
  networks, refrigeration cycles and renewable energy systems.

### Typical applications

1. **Heat exchanger design and optimisation** — double-pipe units,
   evaporators and condensers, with pinch-point analysis to trade off
   thermal effectiveness against surface area.
2. **Pumps, turbines and compressors** — work output, isentropic efficiency
   and energy transfer for power generation, refrigeration and HVAC.
3. **Energy efficiency analysis** — exergy destruction accounting to locate
   where useful work is actually being lost.
4. **Thermodynamic cycle simulation** — Rankine, regenerative Rankine,
   Brayton, ORC and refrigeration cycles.

---

## Resources

- [**PyPI package**](https://pypi.org/project/ThermoSim/)
- [**Source code**](https://github.com/Nouman090/ThermoSim)
- [**Report an issue**](https://github.com/Nouman090/ThermoSim/issues)
- [**Wiki**](https://github.com/Nouman090/ThermoSim/wiki)
- [**Changelog**](https://github.com/Nouman090/ThermoSim/blob/main/CHANGELOG.md)

---

## Contributing

Contributions are welcome. Fork the repository and open a pull request; for
major changes, please open an issue to discuss it first.

---

## License

MIT — see [LICENSE](https://github.com/Nouman090/ThermoSim/blob/main/LICENSE).

---

## Acknowledgements

Created and maintained by
[Md. Waheduzzaman Nouman](https://github.com/Nouman090), for educational and
research use.
