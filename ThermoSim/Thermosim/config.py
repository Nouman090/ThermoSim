"""
config.py
---------
Stores the dead-state (reference environment) used for exergy calculations.
Every other file reads from here instead of from a class variable.
"""

dead_states = {'T0': 298.15, 'P0': 101325}


def set_dead_state(T0=298.15, P0=101325):
    """Update the global dead-state values."""
    global dead_states
    dead_states = {'T0': T0, 'P0': P0}