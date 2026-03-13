
# ThermoSim

![PyPI](https://github.com/Nouman090/ThermoSim/blob/main/docs/ThermoSim%20Logo%202.jpg?raw=true)


**ThermoSim** is a Python package for simulating and analyzing thermodynamic model.
![ThermoSim_logo](https://github.com/Nouman090/ThermoSim/blob/main/docs/ThermoSim%20Logo%202.jpg?raw=true)


## **What is ThermoSim?**

**ThermoSim**, is designed to simulate and analyze various thermodynamic systems and components. Object-Oriented Programming is used in this module. It can model complex systems involving fluids such as water, air, refrigerants like isobutane and so on. The module supports a range of thermodynamic processes, including pumps, turbines, heat exchangers,pipe,expansion valve and other essential components commonly found in energy systems, refrigeration cycles, and heat transfer applications.


## 🔧 Installation

Install the latest release from PyPI:

```bash
pip install thermosim
```

---

## 🚀 Quick Example

<<<<<<< HEAD
```python
import ThermoSim

# Initialize the thermodynamic model
model = ThermoSim.ThermodynamicModel()

# Define fluid state points
model.add_point('water', StatePointName='1', P=6.09e5, T=158+273.15, Mass_flowrate=555.9)
model.add_point('water', StatePointName='2', P=6.09e5, T=None, Mass_flowrate=555.9)

# Add components (e.g., pump, heat exchanger)
pump = model.Pump(model, 'Pump', In_state='1', Out_state='2', n_isen=0.75, Calculate=True)

print(Model)
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
