"""
base_component.py
-----------------
Abstract base class that every component inherits from.
Provides shared helpers so you don't repeat code in every component.
"""

import warnings
from abc import ABC, abstractmethod


class Component(ABC):
    """Base class for all thermodynamic components."""

    def __init__(self, Model, ID, Calculate=False):
        self.Model = Model
        self.ID = ID
        self.Solution_Status = False
        self.Ex_D = "Not Calculated"
        # Register in the model automatically
        self.Model.Component[ID] = self
        if Calculate:
            self.Cal()

    # every child MUST implement Cal() and __str__()
    @abstractmethod
    def Cal(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    # -------------------------------------------------------------- #
    #  shared helpers
    # -------------------------------------------------------------- #
    def _resolve_mass_flowrate(self, pt_in, pt_out):
        """
        Figure out the mass-flow rate from the inlet / outlet points.
        Returns (mass_flowrate, used_default_flag).
        """
        m_in  = pt_in.Mass_flowrate
        m_out = pt_out.Mass_flowrate

        if m_out is None and m_in is not None:
            pt_out.Mass_flowrate = m_in
            return m_in, False

        if m_out is not None and m_in is None:
            pt_in.Mass_flowrate = m_out
            return m_out, False

        if m_out == m_in:
            if m_out is None:
                warnings.warn(
                    f"Mass-flow rate through {self.ID} set to 1 kg/s "
                    f"(none was given)."
                )
                return 1.0, True          # flag = True → used default
            return m_out, False

        if m_out != m_in:
            raise ValueError(
                f"Mass-flow rate mismatch in {self.ID}: "
                f"inlet={m_in}, outlet={m_out}"
            )

        raise ValueError(f"Mass-flow rate of {self.ID} is not given.")

    def _update_model_points(self, **mapping):
        """
        Write component state points back into Model.Point.
        Usage:  self._update_model_points(In_state=self.In,
                                          Out_state=self.Out)
        """
        for attr_name, point_obj in mapping.items():
            state_name = getattr(self, attr_name, None)
            if state_name is not None and point_obj is not None:
                self.Model.Point[state_name] = point_obj