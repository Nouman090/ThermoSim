from setuptools import setup, find_packages

setup(
    name='thermocycle',
    version='3.0.0',
    description='A Python package for thermodynamic cycle modelling and analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Your Name',              # ← change this
    author_email='your@email.com',   # ← change this
    url='https://github.com/YOUR_USERNAME/YOUR_REPO_NAME',  # ← change this
    packages=find_packages(),
    install_requires=[
        'CoolProp',
        'scipy',
        'matplotlib',
        'numpy',
        'pandas',
        'pymoo',
    ],
    extras_require={
        'test': ['pytest'],
    },
    python_requires='>=3.8',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering',
    ],
    package_data={
        'thermocycle': ['T66.json'],
    },
)