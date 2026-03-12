from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()


setup(
      name="ThermoSim",
      version="2.4.0",
      author="Md. Waheduzzaman Basunia Nouman",
      description="A simulation package for thermodynamic systems",
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages= find_packages(),
      url='https://github.com/Nouman090/ThermoSim', 
    project_urls={
        "Bug Tracker": "https://github.com/Nouman090/ThermoSim/issues",
        "Source": "https://github.com/Nouman090/ThermoSim",
        "Wiki": "https://github.com/Nouman090/ThermoSim/wiki"
    },
      license= "MIT",
      install_requires = [
          "numpy",
          "scipy",
          "matplotlib",
          "CoolProp",
          "pandas",
          "CoolProp",
          "pymoo"
          ],
      
      )