"""
model.py
--------
The main ThermodynamicModel class.
This is the "brain" — it holds all state points and components
and provides convenience methods to build and solve cycles.
"""

import warnings
import pandas as pd

from . import config
from .state import Prop


class ThermodynamicModel:
    """Top-level container for a thermodynamic cycle model."""

    # ---------------------------------------------------------- #
    #  constructor  (FIXED: was "init" before → now "__init__")
    # ---------------------------------------------------------- #
    def __init__(self):
        self.Point     = {}      # name → Prop
        self.Component = {}      # name → component object
        self.Connection = []
        self.Loops     = {}

    # ---------------------------------------------------------- #
    #  dead state
    # ---------------------------------------------------------- #
    def set_dead_state(self, T0=298.15, P0=101325):
        config.set_dead_state(T0, P0)

    # ---------------------------------------------------------- #
    #  factory for Prop objects  (so components can call
    #  self.Model.Prop(...)  exactly as before)
    # ---------------------------------------------------------- #
    def Prop(self, fluid, StatePointName='_tmp',
             Mass_flowrate=None, **properties):
        """Create and return a new Prop (state-point) object."""
        return Prop(fluid, StatePointName,
                    Mass_flowrate=Mass_flowrate, **properties)

    # ---------------------------------------------------------- #
    #  add state points
    # ---------------------------------------------------------- #
    def add_point(self, fluid, StatePointName,
                  Mass_flowrate=None, **properties):
        node = Prop(fluid, StatePointName,
                    Mass_flowrate=Mass_flowrate, **properties)
        self.Point[node.StatePointName] = node

    def add_Connection(self, Comp_1, Comp_1_id, Comp_2, Comp_2_id,
                       fluid, StatePointName,
                       Mass_flowrate=None, **properties):
        node = Prop(fluid, StatePointName,
                    Mass_flowrate=Mass_flowrate, **properties)
        self.Point[node.StatePointName] = node
        self.Connection.append([
            Comp_1.ID, Comp_1_id,
            Comp_2.ID, Comp_2_id,
            StatePointName
        ])
        self.Component[Comp_1.ID] = Comp_1
        self.Component[Comp_2.ID] = Comp_2

    # ---------------------------------------------------------- #
    #  iterative solver
    # ---------------------------------------------------------- #
    def Solve(self, max_iter=10, verbose=False):
        """Call Cal() on every component repeatedly until converged."""
        for iteration in range(max_iter):
            # snapshot current enthalpies
            old = {k: v.H for k, v in self.Point.items()
                   if v.H is not None}

            for key, comp in self.Component.items():
                try:
                    comp.Cal()
                except Exception as e:
                    if verbose:
                        print(f"  iter {iteration} | {key} failed: {e}")

            # convergence check
            converged = True
            for k, v in self.Point.items():
                if k in old and v.H is not None:
                    if abs(v.H - old[k]) > 1e-6:
                        converged = False
                        break
            if converged:
                if verbose:
                    print(f"Converged in {iteration + 1} iterations.")
                return True

        if verbose:
            print(f"Did NOT converge after {max_iter} iterations.")
        return False

    # ---------------------------------------------------------- #
    #  pretty-print state points
    # ---------------------------------------------------------- #
    def Point_print(self, header=None):
        default = ['fluid', 'Mass_flowrate', 'StatePointName',
                    'P', 'T', 'H', 'Q', 'ex']
        header = header or default

        rows = [vars(v).copy() for v in self.Point.values()]
        df = pd.DataFrame(rows)

        # unit conversions for display
        if 'T' in df.columns:
            df['T'] = df['T'] - 273.15
        cols = [c for c in header if c in df.columns]
        df = df[cols]
        print(df.to_string(index=False))
        print("\nAvailable columns:", list(pd.DataFrame(rows).columns))
        return df

    # ---------------------------------------------------------- #
    #  pretty-print components
    # ---------------------------------------------------------- #
    def Component_print(self):
        for name, comp in self.Component.items():
            try:
                print(comp, '\n')
            except Exception:
                print(f"{name}: not solved\n")

    # ---------------------------------------------------------- #
    #  model summary  (FIXED: fresh calculation each time)
    # ---------------------------------------------------------- #
    def ModelSummary(self):
        # import here to avoid circular import at module level
        from .turbomachinery   import Turbine, Pump
        from .heat_exchangers  import HeatExchanger

        ex_d_list   = []
        energy_list = []
        Power_Out   = 0
        Power_In    = 0
        Q_in        = 0
        Q_out       = 0
        Total_Ex_d  = 0

        for name, comp in self.Component.items():
            # exergy
            if isinstance(comp.Ex_D, (int, float)):
                Total_Ex_d += comp.Ex_D
            ex_d_list.append([comp.ID, comp.Ex_D])

            # energy
            if isinstance(comp, Turbine):
                Power_Out += comp.work
                energy_list.append([comp.ID + " power out", comp.work])

            elif isinstance(comp, Pump):
                Power_In += comp.work
                energy_list.append([comp.ID + " power in", comp.work])

            elif isinstance(comp, HeatExchanger):
                if comp.HeatAdded is True:
                    Q_in += comp.Q
                    energy_list.append([comp.ID + " heat added", comp.Q])
                elif comp.HeatAdded is False:
                    Q_out += comp.Q
                    energy_list.append([comp.ID + " heat rejected", comp.Q])
                else:
                    energy_list.append([comp.ID + " heat exchanged", comp.Q])

        Net_power  = Power_Out - Power_In
        Efficiency = (Net_power / Q_in * 100) if Q_in != 0 else 0

        energy_list.append(["Net power",        Net_power])
        energy_list.append(["Efficiency (%)",   Efficiency])
        ex_d_list.append(["Total Ex. destr.", Total_Ex_d])

        self.Net_power  = Net_power
        self.Efficiency = Efficiency
        self.Total_Ex_d = Total_Ex_d
        self.Power_Out  = Power_Out
        self.Power_In   = Power_In
        self.Q_in       = Q_in
        self.Q_out      = Q_out

        self.Energy = pd.DataFrame(energy_list,
                                    columns=['Component', 'Energy'])
        self.Ex_D_df = pd.DataFrame(ex_d_list,
                                     columns=['Component', 'Ex_destruction'])

        print("═" * 50)
        print("  ENERGY ANALYSIS")
        print("═" * 50)
        print(self.Energy.to_string(index=False))
        print()
        print("═" * 50)
        print("  EXERGY ANALYSIS")
        print("═" * 50)
        print(self.Ex_D_df.to_string(index=False))
        print()
        # ---------------------------------------------------------- #
    #  save / load
    # ---------------------------------------------------------- #
    def save_model(self, filename):
        """
        Save the current model state to a JSON file.

        Saves every state point (fluid, mass flow, P, T, H, S, Q)
        and component definitions so the model can be reconstructed.

        Parameters
        ----------
        filename : str
            Path to the output JSON file, e.g. 'my_cycle.json'
        """
        import json

        data = {
            'dead_state': config.dead_states,
            'points': {},
            'components': {},
        }

        # --- save all state points -----------------------------------
        for name, pt in self.Point.items():
            data['points'][name] = {
                'fluid':          pt.fluid,
                'Mass_flowrate':  pt.Mass_flowrate,
                'P': pt.P,
                'T': pt.T,
                'H': pt.H,
                'S': pt.S,
                'Q': pt.Q if isinstance(pt.Q, (int, float, type(None))) else None,
            }

        # --- save component definitions ------------------------------
        for name, comp in self.Component.items():
            comp_data = {'type': type(comp).__name__}
            # grab common attributes
            for attr in ['In_state', 'Out_state', 'In_states',
                         'Out_states', 'split_fractions',
                         'n_isen', 'n_mech', 'Compressibility',
                         'Hot_In_state', 'Hot_Out_state',
                         'Cold_In_state', 'Cold_Out_state',
                         'PPT', 'HEX_type', 'HeatAdded',
                         'effectiveness', 'Q', 'UA',
                         'Pressure_drop', 'Temperature_drop',
                         'Out_vap_state', 'Out_liq_state',
                         'T_melt', 'per_loss', 'Charge',
                         'Charging_time', 'Discharging_time',
                         'Capacity', 'pressure_ratio']:
                if hasattr(comp, attr):
                    val = getattr(comp, attr)
                    # only save JSON-serialisable values
                    if isinstance(val, (int, float, str, bool,
                                        list, type(None))):
                        comp_data[attr] = val
            data['components'][name] = comp_data

        with open(filename, 'w') as f:
            json.dump(data, f, indent=2, default=str)

        print(f"Model saved to '{filename}'")

    def load_model(self, filename):
        """
        Load state points from a JSON file previously created by
        save_model().  Points are restored so you can re-create
        components and solve.

        Parameters
        ----------
        filename : str
            Path to the JSON file.
        """
        import json

        with open(filename, 'r') as f:
            data = json.load(f)

        # restore dead state
        ds = data.get('dead_state', {})
        self.set_dead_state(
            T0=ds.get('T0', 298.15),
            P0=ds.get('P0', 101325),
        )

        # restore points
        for name, pt_data in data.get('points', {}).items():
            fluid = pt_data['fluid']
            mf    = pt_data.get('Mass_flowrate')

            # pick the best two known properties
            candidates = {}
            for key in ['P', 'T', 'H', 'S', 'Q']:
                val = pt_data.get(key)
                if val is not None and isinstance(val, (int, float)):
                    candidates[key] = val

            # prefer P + one other
            chosen = {}
            if 'P' in candidates:
                chosen['P'] = candidates.pop('P')
            # pick the next available
            for key in ['H', 'T', 'S', 'Q']:
                if key in candidates and len(chosen) < 2:
                    chosen[key] = candidates[key]

            if len(chosen) == 2:
                self.add_point(fluid, name, Mass_flowrate=mf, **chosen)
            else:
                # if we can't get two props, just store what we have
                self.add_point(fluid, name, Mass_flowrate=mf, **chosen)
                warnings.warn(
                    f"Point '{name}': only {len(chosen)} properties "
                    f"could be restored."
                )

        print(f"Model loaded from '{filename}'")
        print(f"  {len(self.Point)} state points restored.")
        print(f"  Component definitions in file: "
              f"{list(data.get('components', {}).keys())}")
        print("  NOTE: Components must be re-created from the script "
              "or re-attached manually.")
    # ---------------------------------------------------------- #
    def __str__(self):
        print("\n" + "=" * 60)
        print("           MODEL SUMMARY")
        print("=" * 60 + "\n")
        self.ModelSummary()
        print("\n--- Component Details ---\n")
        self.Component_print()
        print("\n--- State Point Table ---\n")
        self.Point_print()
        return ''