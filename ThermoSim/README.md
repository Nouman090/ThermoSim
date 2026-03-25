
# ThermoSim

![PyPI](https://github.com/Nouman090/ThermoSim/blob/main/docs/ThermoSim%20Logo%202.jpg?raw=true)


**ThermoSim** is a Python package for simulating and analyzing thermodynamic model.
![ThermoSim_logo](https://github.com/Nouman090/ThermoSim/blob/main/docs/ThermoSim%20Logo%202.jpg?raw=true)

# ThermoCycle

A Python package for modelling and analysing thermodynamic power cycles.

## Features

- **State point calculation** using CoolProp (water, refrigerants, air, etc.)
- **Components**: Turbine, Pump, Compressor, Heat Exchanger, Mixer, Splitter, 
  Expansion Valve, Separator, Pipe, TES, Source, Sink
- **Exergy analysis** with dead-state reference
- **Cycle plots**: T-s, P-h, h-s diagrams with saturation domes
- **Pinch analysis** for heat exchangers
- **Sensitivity analysis**: single-parameter, multi-output, and 2D contour sweeps
- **Save/Load** model state to JSON
- **48+ unit tests** with pytest

## Installation

```bash
pip install CoolProp scipy matplotlib numpy pandas pymoo pytest
```
---

## 🚀 Quick Example

<<<<<<< HEAD
```python
from thermocycle import (
    ThermodynamicModel, Turbine, Pump, HeatExchanger,
)
from thermocycle.plotting import CyclePlotter

Model = ThermodynamicModel()
Model.set_dead_state()

Model.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
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

##Test Run
```python
pytest tests/ -v

```
---

## 📚 Resources

- 🧾 [**PyPI Package**](https://pypi.org/project/ThermoSim/)
- 🛠 [**Source Code**](https://github.com/Nouman090/ThermoSim)
- ❓ [**Report Issues**](https://github.com/Nouman090/ThermoSim/issues)
- 📘 [**Wiki**](https://github.com/Nouman090/ThermoSim/wiki)

---

## ✨ Features

- Support for different working and heating fluids
- Handles mass and energy balances automatically
- Designed for academic and research-grade modeling

---

## 🤝 Contributing

You're welcome to contribute! Fork the repository and open a pull request. For major changes, please discuss via an issue first.

---

## 📄 License

This project is licensed under the MIT License.  
See the [LICENSE](https://github.com/Nouman090/ThermoSim/blob/main/LICENSE) file for more details.

---

## 🙌 Acknowledgements

Created and maintained by [Md. Waheduzzaman Nouman](https://github.com/Nouman090), for educational and research use.
=======
This documentation is intended for:

* **Engineering Students**: Those studying thermodynamics, energy systems, and heat transfer. The module provides a practical tool to simulate real-world energy systems and understand thermodynamic concepts.

* **Researchers**: Professionals and researchers working in the field of thermodynamics and energy efficiency can use this module for modeling, optimization, and analysis of complex systems.

* **Energy System Designers**: Engineers involved in designing and optimizing thermodynamic systems such as power plants, heat exchangers, refrigeration cycles, and renewable energy systems.

## **Real-World Applications**

1. **Heat Exchanger Design and Optimization**: The module simulates various types of heat exchangers (e.g., double-pipe, evaporator, condenser), helping engineers optimize thermal efficiency and energy usage in industrial applications.

2. **Pumps and Turbines**: It can model pumps and turbines used in power generation, refrigeration, and HVAC systems, providing insights into performance metrics like work output, efficiency, and energy transfer  
     
3. **Energy Efficiency Analysis**: By integrating components like expansion valves and PCM (Phase Change Materials), the model supports the design of energy-efficient systems in heating, cooling, and refrigeration sectors.  
     
4. **Simulation of Thermodynamic Cycles**: The module supports the simulation of thermodynamic cycles, including Rankine and refrigeration cycles, helping in the evaluation of system performance, energy conservation, and operational optimization.

# License

MIT
>>>>>>> c93e50c82d03eece8612d6db09b20b1bbd2ce23c
