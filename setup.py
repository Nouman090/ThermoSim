from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()


setup(
      name="ThermoSim",
      version="2.3.2",
      author="Md. Waheduzzaman Basunia Nouman",
      description="A simulation package for thermodynamic systems",
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages= find_packages(),
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