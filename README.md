# ThermoSim

![PyPI](https://github.com/Nouman090/ThermoSim/blob/main/docs/ThermoSim%20Logo%202.jpg?raw=true)


**ThermoSim** is a Python package for simulating and analyzing thermodynamic model.

## 🔧 Installation

Install the latest release from PyPI:

```bash
pip install thermosim
```

---

## 🚀 Quick Example

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
