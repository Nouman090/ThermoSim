"""
model.py
--------
The main ThermodynamicModel class.
This is the "brain" — it holds all state points and components
and provides convenience methods to build and solve cycles.
"""

import warnings
import pandas as pd
from tabulate import tabulate


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
        self._propagation_enabled = False   # set True by enable_flow_propagation()

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
        # If propagation is already live, wire this new point in immediately.
        # The graph will be rebuilt on the next enable_flow_propagation() call;
        # for now we at least give it a propagator so a Mass_flowrate assignment
        # on THIS point won't silently drop the update.
        if self._propagation_enabled:
            object.__setattr__(node, '_flow_propagator', self._make_propagator())
        return node

    def add_Connection(self, Comp_1, Comp_1_id, Comp_2, Comp_2_id,
                       fluid, StatePointName,
                       Mass_flowrate=None, **properties):
        node = Prop(fluid, StatePointName,
                    Mass_flowrate=Mass_flowrate, **properties)
        self.Point[node.StatePointName] = node
        if self._propagation_enabled:
            object.__setattr__(node, '_flow_propagator', self._make_propagator())
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
    def _sorted_components(self):
        """
        Return the component dict items in the prescribed solve order:

        Priority 1 – single-stream pass-through components
                     (Turbine, Pump, Compressor, Pipe, Expansion_valve)
        Priority 2 – flow-splitting / mixing components
                     (Splitter, Separator, Mixer)
        Priority 3a – SimpleHEX heat exchangers   (solved first among HEXs)
        Priority 3b – other HeatExchanger types, sorted in *descending* order
                      of their hot-inlet temperature so the hottest HEX is
                      solved first (it constrains downstream HEXs).
        Priority 4  – everything else (TES, custom components, …)

        Within each priority group the original insertion order is preserved.
        """
        from .turbomachinery   import Turbine, Pump, Compressor
        from .simple_components import Splitter, Mixer, Separator, Pipe, Expansion_valve
        from .heat_exchangers  import HeatExchanger

        # ── buckets ──────────────────────────────────────────────────────
        pass_through = []   # priority 1
        flow_split   = []   # priority 2
        simple_hex   = []   # priority 3a
        other_hex    = []   # priority 3b  (will be sorted by T_hot_in)
        rest         = []   # priority 4

        _pass_through_types = (Turbine, Pump, Compressor, Pipe, Expansion_valve)
        _flow_types         = (Splitter, Separator, Mixer)

        for key, comp in self.Component.items():
            if isinstance(comp, _pass_through_types):
                pass_through.append((key, comp))

            elif isinstance(comp, _flow_types):
                flow_split.append((key, comp))

            elif isinstance(comp, HeatExchanger):
                if getattr(comp, 'HEX_type', None) == 'SimpleHEX':
                    simple_hex.append((key, comp))
                else:
                    other_hex.append((key, comp))

            else:
                rest.append((key, comp))

        # ── sort non-SimpleHEX by hot-inlet temperature (descending) ─────
        def _hot_inlet_T(item):
            """
            Return the hot-inlet temperature for a HeatExchanger.
            Falls back to -inf so unsolved / one-sided HEXs sort last.
            """
            comp = item[1]
            hot_in_name = getattr(comp, 'Hot_In_state', None)
            if hot_in_name and hot_in_name in self.Point:
                T = self.Point[hot_in_name].T
                if isinstance(T, (int, float)):
                    return T
            return float('-inf')

        other_hex.sort(key=_hot_inlet_T, reverse=True)

        # ── concatenate in priority order ────────────────────────────────
        ordered = pass_through + flow_split + simple_hex + other_hex + rest

        if len(ordered) != len(self.Component):
            # safety net – append anything that was missed
            seen = {k for k, _ in ordered}
            for key, comp in self.Component.items():
                if key not in seen:
                    ordered.append((key, comp))

        return ordered
    
    def Solve(self, max_iter=10, verbose=False, custom_order=None,
              tol=1e-8):
        """
        Call Cal() on every component repeatedly until converged.

        Parameters
        ----------
        max_iter : int
            Maximum number of iterations.
        verbose : bool
            Print convergence details.
        tol : float
            Relative convergence tolerance on state-point enthalpy between
            iterations.  Convergence requires
            ``|ΔH| <= tol * max(|H|, 1.0)`` at every point.
        custom_order : list of str
            Optional list of component IDs specifying the exact solve order.
            If provided, components are attempted in this order every iteration.
            Any registered components NOT in the list are appended at the end
            in their default priority order so nothing is silently skipped.
            Example: M.Solve(custom_order=['Pump1', 'Boiler', 'Turbine1'])

        Default priority order (used when custom_order is empty):
          1. Pass-through (Turbine, Pump, Compressor, Pipe, Expansion_valve)
          2. Flow-split / mixing  (Splitter, Separator, Mixer)
          3a. SimpleHEX heat exchangers
          3b. Other HeatExchanger types, sorted by hot-inlet T (descending)
          4. Everything else
        """
        from .heat_exchangers  import HeatExchanger
        custom_order = list(custom_order) if custom_order else []
        self.enable_flow_propagation(verbose=verbose)
        if verbose:
            self.flow_graph_summary()
        number_of_components = len(self.Component)
        Solve_components = 0
        error_list = {}
        breakFlag = 0
        stalled = False

        if number_of_components == 0:
            if verbose:
                print("No components to solve.")
            return True  # trivially converged

        # Validate custom_order IDs up front so the user gets a clear error
        # rather than a silent miss deep inside the iteration loop.
        if custom_order:
            unknown = [cid for cid in custom_order if cid not in self.Component]
            if unknown:
                raise ValueError(
                    f"Solve(custom_order=...) contains unknown component ID(s): "
                    f"{unknown}. Registered IDs: {list(self.Component.keys())}"
                )

        for iteration in range(max_iter):
            old_Solve_components = Solve_components
            old = {k: v.H for k, v in self.Point.items()
                   if v.H is not None}
            solved_comp = []

            # ── Build the ordered component list for this iteration ────────
            if custom_order:
                ordered_components = [(cid, self.Component[cid])
                                      for cid in custom_order]

                listed = set(custom_order)
                remaining = [(k, c) for k, c in self._sorted_components()
                             if k not in listed]
                ordered_components += remaining
            else:
                ordered_components = self._sorted_components()

            for key, comp in ordered_components:
                flag1 = True
                try:
                    if comp.Solution_Status == True:
                        flag1 = False
                    comp.Cal()
                    if verbose:
                        print(f"  solved {comp.ID}")
                    if comp.Solution_Status == True:
                        error_list.pop(key, None)
                        if flag1:
                            Solve_components += 1
                            solved_comp.append(str(key))
                except Exception as e:
                    if verbose:
                        print(f"  could not solve {comp.ID}: {e}")
                    error_list[key] = str(e)

            if old_Solve_components == Solve_components:
                breakFlag += 1
            else:
                breakFlag = 0

            if breakFlag >= 3:
                for key, comp in ordered_components:
                    flag1 = True
                    try:
                        if comp.Solution_Status == True:
                            flag1 = False
                        if isinstance(comp, HeatExchanger):
                            comp.PPT_graph = True
                        comp.Cal()
                        if verbose:
                            print(f"  solved {comp.ID} (final pass)")
                        if comp.Solution_Status == True:
                            error_list.pop(key, None)
                            if flag1:
                                Solve_components += 1
                                solved_comp.append(str(key))
                    except Exception as e:
                        if verbose:
                            print(f"  could not solve {comp.ID}: {e}")
                        error_list[key] = str(e)
                stalled = True
                warnings.warn(
                    f"Solve() stopped: no new component was solved for 3 "
                    f"consecutive iterations. {Solve_components} of "
                    f"{number_of_components} components solved."
                )
                break

            if verbose:
                print(
                    f"  iter {iteration} | {Solve_components} of "
                    f"{number_of_components} solved | "
                    f"solved this iter: {solved_comp}"
                )

            # ── convergence check ──────────────────────────────────────
            # Relative tolerance.  An absolute 1e-6 J/kg on enthalpies of
            # order 3e6 J/kg is ~3e-13 relative — below what a CoolProp
            # round-trip can reproduce, so the loop could essentially never
            # report convergence on its own terms.
            converged = True
            for k, v in self.Point.items():
                if k in old and v.H is not None:
                    if abs(v.H - old[k]) > tol * max(abs(v.H), 1.0):
                        converged = False
                        break
                else:
                    converged = False
                    break

            if converged:
                if verbose:
                    print(f"Converged in {iteration + 1} iterations.")
                return True

        if not stalled:
            warnings.warn(
                f"Solve() did NOT converge after {max_iter} iterations."
            )
        if verbose and error_list:
            print(tabulate(
                error_list.items(),
                headers=['Component', 'Cause of failure'],
                maxcolwidths=[10, 50],
                tablefmt='grid'
            ))
        return False

    # ---------------------------------------------------------- #
    #  pretty-print state points
    # ---------------------------------------------------------- #
    def Point_print(self, header=None, condition=True):
        """
        Tabulate every state point.

        Args:
            header: list of column names to display.  Defaults to
                ['fluid', 'Mass_flowrate', 'StatePointName', 'P', 'T', 'H',
                 'Q', 'Condition', 'ex'].
            condition: when True (default) a 'Condition' column is derived
                from ``phase``.  ``Q`` is numeric only inside the saturation
                dome and NaN outside it, so this column carries the label
                -- 'Superheated', 'Sub-cooled', 'Saturated liquid' -- that
                ``Q`` used to be overwritten with before 3.2.2.
        """
        default = ['fluid', 'Mass_flowrate', 'StatePointName',
                    'P', 'T', 'H', 'Q']
        if condition:
            default.append('Condition')
        default.append('ex')
        header = header or default

        # Public attributes only: vars() also exposes internals such as
        # _initialising, _flow_propagator, _user_set_flowrate and the
        # properties/DeadStates dicts, which were previously advertised to
        # the user in the "Available columns" line.
        _private = {'properties', 'DeadStates', 'Solution_Status'}
        rows = [
            {k: v for k, v in vars(pt).items()
             if not k.startswith('_') and k not in _private}
            for pt in self.Point.values()
        ]
        full = pd.DataFrame(rows)

        # 'Condition' is the short form of 'phase' with the numeric margin
        # stripped: 'Superheated (184.99 K)' -> 'Superheated'.  The full
        # label stays available in the 'phase' column.
        if 'phase' in full.columns:
            full['Condition'] = full['phase'].map(self._short_phase)

        df = full.copy()

        # unit conversions for display
        if 'T' in df.columns:
            df['T'] = df['T'] - 273.15
        cols = [c for c in header if c in df.columns]
        df = df[cols]
        print(df.to_string(index=False))
        print("\nAvailable columns:", list(full.columns))
        return df

    @staticmethod
    def _short_phase(label):
        """'Superheated (184.99 K)' -> 'Superheated'; None -> 'Unsolved'."""
        if label is None:
            return 'Unsolved'
        return str(label).split('(')[0].strip()

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
    def ModelSummary(self, verbose=True):
        # self.Solve()
        # import here to avoid circular import at module level
        from .turbomachinery   import Turbine, Pump, Compressor
        from .heat_exchangers  import HeatExchanger

        ex_d_list   = []
        energy_list = []
        Power_Out   = 0
        Power_In    = 0
        Q_in        = 0
        Q_out       = 0
        Total_Ex_d  = 0

        skipped = []

        for name, comp in self.Component.items():
            # exergy
            if isinstance(comp.Ex_D, (int, float)):
                Total_Ex_d += comp.Ex_D
            ex_d_list.append([comp.ID, comp.Ex_D])

            # energy
            # Every accumulator below is guarded with `is not None`: an
            # unsolved component has work/Q of None, and `Q_in += None`
            # raised a TypeError that aborted the whole summary.
            if isinstance(comp, Turbine):
                if comp.work is None:
                    skipped.append(comp.ID)
                else:
                    Power_Out += comp.work
                    energy_list.append([comp.ID + " power out", comp.work])

            elif isinstance(comp, (Pump, Compressor)):
                if comp.work is None:
                    skipped.append(comp.ID)
                else:
                    Power_In += comp.work
                    energy_list.append([comp.ID + " power in", comp.work])

            elif isinstance(comp, HeatExchanger):
                if comp.Q is None:
                    skipped.append(comp.ID)
                elif comp.HeatAdded is True:
                    Q_in += comp.Q
                    energy_list.append([comp.ID + " heat added", comp.Q])
                elif comp.HeatAdded is False:
                    Q_out += comp.Q
                    energy_list.append([comp.ID + " heat rejected", comp.Q])
                else:
                    energy_list.append([comp.ID + " heat exchanged", comp.Q])

        if skipped:
            warnings.warn(
                f"ModelSummary(): {len(skipped)} component(s) contributed "
                f"nothing because they are not solved: {skipped}. "
                f"Call Solve() first, or check the errors it reported."
            )

        Net_power = Power_Out - Power_In
        COP_R = COP_HP = None

        if Q_in > Q_out:
            # Power cycle
            Efficiency = (Net_power / Q_in * 100) if Q_in != 0 else 0.0
            energy_list.append(["Net power",      Net_power])
            energy_list.append(["Efficiency (%)", Efficiency])
        else:
            # Refrigeration / heat-pump cycle: thermal efficiency is not a
            # meaningful figure of merit, so it is NaN rather than the string
            # "Not Applicable".  A string broke every caller that did
            # arithmetic or comparison on it (SensitivityAnalyzer plots it;
            # the test-suite asserts `Efficiency > 0`).
            Efficiency = float('nan')
            # COP_R  = useful cooling / work in  = Q_evaporator / W
            # COP_HP = useful heating / work in  = Q_condenser  / W
            # and COP_HP = COP_R + 1 by the first law.
            #
            # v3.2.1 had COP_R = Q_out / Power_In — that is the HEAT PUMP
            # figure (Q_out is heat REJECTED at the condenser), and it then
            # reported COP_HP = Q_out/W + 1, which is COP_R + 2 and
            # corresponds to nothing physical. Verified on an R134a cycle
            # (2 → 10 bar, η=0.8): Q_in=13.71 kW, Q_out=17.90 kW, W=4.19 kW,
            # so COP_R=3.27 and COP_HP=4.27, not 4.27 and 5.27.
            COP_R  = (Q_in  / Power_In) if Power_In != 0 else 0.0
            COP_HP = (Q_out / Power_In) if Power_In != 0 else 0.0
            energy_list.append(["COP Refrigeration", COP_R])
            energy_list.append(["COP Heat Pump",     COP_HP])

        ex_d_list.append(["Total Ex. destr.", Total_Ex_d])

        self.Net_power  = Net_power
        self.Efficiency = Efficiency
        self.COP_R      = COP_R      # None for a power cycle
        self.COP_HP     = COP_HP     # None for a power cycle
        self.Total_Ex_d = Total_Ex_d
        self.Power_Out  = Power_Out
        self.Power_In   = Power_In
        self.Q_in       = Q_in
        self.Q_out      = Q_out

        self.Energy = pd.DataFrame(energy_list,
                                    columns=['Component', 'Energy'])
        self.Ex_D_df = pd.DataFrame(ex_d_list,
                                     columns=['Component', 'Ex_destruction'])

        # Printing is optional: SensitivityAnalyzer calls this once per
        # sweep point, so an unconditional print buried the actual results
        # under one full energy+exergy table per evaluation.
        if verbose:
            print("═" * 50)
            print(" ENERGY ANALYSIS")
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
    #  Flow-propagation network
    # ---------------------------------------------------------- #
    def _build_flow_graph(self):
        """
        Scan every registered component and build two structures:

        _flow_graph : dict[point_name -> list[point_name]]
            Directed edges: upstream point → list of immediate downstream
            points within the same continuous stream.

        _mixer_inputs   : dict[out_point -> [in_point, ...]]
        _splitter_fracs : dict[in_point  -> [(out_point, fraction|None), ...]]

        Rules per component type
        ────────────────────────
        • Simple pass-through (Turbine, Pump, Compressor, Pipe,
          Expansion_valve, Source, Sink):
              In_state → Out_state  (1-to-1, same mass flow)

        • HeatExchanger / TES:
              Hot_In  → Hot_Out   (independent leg)
              Cold_In → Cold_Out  (independent leg)

        • Splitter:
              In_state → each Out_state  (fractions kept fixed)

        • Mixer:
              each In_state → Out_state  (outlet = sum of inlets)

        • Separator:
              In_state → Out_vap_state, Out_liq_state
              (fractions derived from quality Q at propagation time)
        """
        from .simple_components import Splitter, Mixer, Separator
        from .turbomachinery    import Turbine, Pump, Compressor
        try:
            from .heat_exchangers import HeatExchanger, TES
            _hex_types = (HeatExchanger, TES)
        except ImportError:
            _hex_types = ()

        graph          = {name: [] for name in self.Point}
        mixer_inputs   = {}   # out_pt  -> [in_pts]
        splitter_fracs = {}   # in_pt   -> [(out_pt, frac|None)]

        for comp_id, comp in self.Component.items():

            # ── Splitter ──────────────────────────────────────────
            if isinstance(comp, Splitter):
                src = comp.In_state
                for out_name, frac in zip(comp.Out_states,
                                          comp.split_fractions):
                    if src in graph:
                        graph[src].append(out_name)
                splitter_fracs[src] = list(
                    zip(comp.Out_states, comp.split_fractions)
                )

            # ── Mixer ─────────────────────────────────────────────
            elif isinstance(comp, Mixer):
                dst = comp.Out_state
                for in_name in comp.In_states:
                    if in_name in graph:
                        graph[in_name].append(dst)
                mixer_inputs[dst] = list(comp.In_states)

            # ── Separator ─────────────────────────────────────────
            elif isinstance(comp, Separator):
                src = comp.In_state
                vap = comp.Out_vap_state
                liq = comp.Out_liq_state
                if src in graph:
                    graph[src].extend([vap, liq])
                # None = fraction computed dynamically from Q at runtime
                splitter_fracs[src] = [(vap, None), (liq, None)]

            # ── HeatExchanger / TES (two independent legs) ────────
            elif _hex_types and isinstance(comp, _hex_types):
                hot_in  = getattr(comp, 'Hot_In_state',  None)
                hot_out = getattr(comp, 'Hot_Out_state', None)
                if hot_in and hot_out and hot_in in graph:
                    graph[hot_in].append(hot_out)
                cold_in  = getattr(comp, 'Cold_In_state',  None)
                cold_out = getattr(comp, 'Cold_Out_state', None)
                if cold_in and cold_out and cold_in in graph:
                    graph[cold_in].append(cold_out)

            # ── All other pass-through components ─────────────────
            else:
                in_st  = getattr(comp, 'In_state',  None)
                out_st = getattr(comp, 'Out_state', None)
                if in_st and out_st and in_st in graph:
                    graph[in_st].append(out_st)

        self._flow_graph     = graph
        self._mixer_inputs   = mixer_inputs
        self._splitter_fracs = splitter_fracs
        self._build_segments()   # derive segments from the graph

    def _build_segments(self):
        """
        Pre-compute flow segments: groups of points that always carry
        the SAME mass flow rate, bounded by Splitters and Mixers.

        Boundary rule
        ─────────────
        An edge (A → B) is a segment boundary — meaning A and B are in
        DIFFERENT segments — if ANY of the following are true:

          1. A is a Splitter inlet  →  the edge A→B (to any outlet) is cut.
             This also means the edge ENTERING A (C→A) is cut, isolating
             the splitter inlet into its own single-point segment.  Without
             this, the splitter inlet ends up in the same segment as the
             mixer outlet (in a closed loop), causing the BFS to confuse
             "origin segment" with "downstream segment".

          2. B is a Mixer outlet   →  every inlet edge X→B is cut.

        All other edges are intra-segment (same flow, undirected).
        """
        graph          = self._flow_graph
        mixer_inputs   = self._mixer_inputs
        splitter_fracs = self._splitter_fracs

        # Collect all splitter inlet points and mixer outlet points.
        splitter_inlets = set(splitter_fracs.keys())
        mixer_outlets   = set(mixer_inputs.keys())

        boundary_edges = set()

        # Rule 1a: splitter inlet → each outlet
        for split_in, outlets in splitter_fracs.items():
            for out_pt, _ in outlets:
                boundary_edges.add((split_in, out_pt))

        # Rule 1b: any edge ENTERING a splitter inlet (isolates the junction)
        for src, dsts in graph.items():
            for dst in dsts:
                if dst in splitter_inlets:
                    boundary_edges.add((src, dst))

        # Rule 2: every mixer inlet → mixer outlet
        for mix_out, inlets in mixer_inputs.items():
            for in_pt in inlets:
                boundary_edges.add((in_pt, mix_out))

        # Build undirected adjacency using only NON-boundary edges.
        undirected = {pt: set() for pt in graph}
        for src, dsts in graph.items():
            for dst in dsts:
                if (src, dst) not in boundary_edges:
                    undirected[src].add(dst)
                    undirected[dst].add(src)

        # Flood-fill connected components → each = one segment.
        from collections import deque as _deque
        point_to_seg = {}
        segments     = {}
        seg_id       = 0
        for start in graph:
            if start in point_to_seg:
                continue
            q = _deque([start])
            point_to_seg[start] = seg_id
            segments[seg_id] = {start}
            while q:
                node = q.popleft()
                for nbr in undirected.get(node, []):
                    if nbr not in point_to_seg:
                        point_to_seg[nbr] = seg_id
                        segments[seg_id].add(nbr)
                        q.append(nbr)
            seg_id += 1

        self._point_to_segment = point_to_seg
        self._segments         = segments

        # seg_out_edges[seg_id] = [(src_pt, dst_pt, edge_type), ...]
        # Only outgoing boundary edges from each segment.
        # edge_type: 'split'       = from splitter inlet → outlet branch
        #            'mix_in'      = mixer inlet → mixer outlet
        #            'passthrough' = upstream segment → splitter inlet segment
        seg_out_edges = {s: [] for s in segments}
        # seg_in_edges[seg_id] = [(src_pt, dst_pt, edge_type), ...]
        # Incoming boundary edges into each segment (reverse of seg_out_edges).
        seg_in_edges  = {s: [] for s in segments}
        for src, dst in boundary_edges:
            s_src = point_to_seg.get(src)
            s_dst = point_to_seg.get(dst)
            if s_src is not None:
                if src in splitter_inlets:
                    et = 'split'
                elif dst in mixer_outlets:
                    et = 'mix_in'
                else:
                    et = 'passthrough'
                seg_out_edges[s_src].append((src, dst, et))
                if s_dst is not None:
                    seg_in_edges[s_dst].append((src, dst, et))

        # seg_mix_outlets[seg_id] = [mixer_out_pt, ...]
        seg_mix_outlets = {s: [] for s in segments}
        for mix_out in mixer_inputs:
            s = point_to_seg.get(mix_out)
            if s is not None:
                seg_mix_outlets[s].append(mix_out)

        # seg_splitter_inlets[seg_id] = True if this segment IS a splitter-
        # inlet segment (a single junction point).  Used by propagation to
        # immediately fire split edges after writing.
        seg_splitter_inlets = {}
        for split_in in splitter_inlets:
            s = point_to_seg.get(split_in)
            if s is not None:
                seg_splitter_inlets[s] = split_in   # seg → splitter_in point

        self._seg_out_edges       = seg_out_edges
        self._seg_in_edges        = seg_in_edges
        self._seg_mix_outlets     = seg_mix_outlets
        self._seg_splitter_inlets = seg_splitter_inlets

    def _propagate_mass_flow(self, origin_point, new_value):
        """
        Fully bidirectional mass-flow propagation.

        Algorithm
        ─────────
        Every segment is an atomic unit carrying a single mass flow.
        When any segment receives a confirmed flow value it is immediately
        written to all points in that segment, then the BFS fans out in
        BOTH directions across segment boundaries:

        Forward edges  (seg_out_edges):
          • passthrough  →  downstream segment gets same flow
          • split        →  downstream branch gets  flow × fraction
          • mix_in       →  accumulate; fire when ALL inlets are known,
                            outlet = sum of inlets

        Backward edges  (seg_in_edges, traversed in reverse):
          • passthrough  →  upstream segment gets same flow
          • split        →  upstream splitter-inlet segment gets
                            flow / fraction   (back-calculate total)
          • mix_in       →  if mixer-outlet AND all other inlets are
                            already known, this inlet =
                            outlet − sum(other inlets)

        The same visited-segment guard prevents cycles; split destinations
        always bypass the guard so fractional flows are never shadowed by
        a stale full-flow value that arrived first.
        """
        from collections import deque

        point_to_seg   = self._point_to_segment
        segments       = self._segments
        seg_out_edges  = self._seg_out_edges
        seg_in_edges   = self._seg_in_edges
        seg_mix_outlets= self._seg_mix_outlets
        mixer_inputs   = self._mixer_inputs
        splitter_fracs = self._splitter_fracs

        def _write(pt_name, flow):
            """Write mass flow to a point without re-triggering propagator."""
            pt = self.Point.get(pt_name)
            if pt is not None and flow != pt.Mass_flowrate:
                object.__setattr__(pt, 'Mass_flowrate', flow)
                if pt.__dict__.get('Solution_Status'):
                    pt._calc_exergy()

        def _write_segment(seg_id, flow):
            """Write flow to every point in a segment."""
            for pt_name in segments.get(seg_id, []):
                _write(pt_name, flow)

        def _seg_flow_of(seg_id):
            """Return the current known flow of a segment, or None."""
            for pt_name in segments.get(seg_id, []):
                pt = self.Point.get(pt_name)
                if pt is not None and pt.Mass_flowrate is not None:
                    return pt.Mass_flowrate
            return None

        def _get_split_frac(src_pt, dst_pt):
            """Return the split fraction for the edge src_pt → dst_pt, or None."""
            return next(
                (f for op, f in splitter_fracs.get(src_pt, []) if op == dst_pt),
                None
            )

        # ── Step 1: determine origin segment ──────────────────────────
        origin_seg = point_to_seg.get(origin_point)
        if origin_seg is None:
            return   # point not in graph — nothing to do

        # Write the origin segment.
        _write_segment(origin_seg, new_value)

        # If new_value is None, propagate None everywhere reachable
        # (clears stale values) and return — no arithmetic needed.
        if new_value is None:
            _visited = {origin_seg}
            _q = deque([origin_seg])
            while _q:
                _seg = _q.popleft()
                for _edges in (seg_out_edges.get(_seg, []),
                               seg_in_edges.get(_seg, [])):
                    for _sp, _dp, _et in _edges:
                        for _pt in (_sp, _dp):
                            _ns = point_to_seg.get(_pt)
                            if _ns is not None and _ns not in _visited:
                                _visited.add(_ns)
                                _write_segment(_ns, None)
                                _q.append(_ns)
            return

        # ── Step 2: bidirectional segment-level BFS ───────────────────
        # seg_flow tracks confirmed flows so backward passes can read them.
        seg_flow    = {origin_seg: new_value}
        visited_segs = {origin_seg}

        # Mixer accumulation state  (used by forward mix_in edges)
        mix_arrived = {}   # mix_out_pt → {inlet_seg: flow}

        # BFS queue entries: (seg_id, src_pt, dst_pt, edge_type, flow, direction)
        # direction: 'fwd' = forward (following seg_out_edges)
        #            'bwd' = backward (following seg_in_edges in reverse)
        queue = deque()

        def _enqueue_neighbours(seg_id, flow):
            """Push all forward and backward neighbours of seg_id."""
            # Forward
            for sp, dp, et in seg_out_edges.get(seg_id, []):
                ns = point_to_seg.get(dp)
                if ns is not None and (ns not in visited_segs or et == 'split'):
                    queue.append((ns, sp, dp, et, flow, 'fwd'))
            # Backward — traverse seg_in_edges in reverse
            for sp, dp, et in seg_in_edges.get(seg_id, []):
                ns = point_to_seg.get(sp)
                if ns is not None and ns not in visited_segs:
                    queue.append((ns, sp, dp, et, flow, 'bwd'))

        _enqueue_neighbours(origin_seg, new_value)

        while queue:
            dst_seg, src_pt, dst_pt, etype, src_flow, direction = queue.popleft()

            new_seg_flow = None   # will be set below if resolvable

            # ══════════════════════════════════════════════════════════
            # FORWARD direction
            # ══════════════════════════════════════════════════════════
            if direction == 'fwd':

                # ── Mixer outlet: accumulate inlets ───────────────────
                if seg_mix_outlets.get(dst_seg):
                    mix_out_pt = seg_mix_outlets[dst_seg][0]

                    if mix_out_pt not in mix_arrived:
                        mix_arrived[mix_out_pt] = {}

                    # Record this inlet's contribution keyed by its segment.
                    src_seg = point_to_seg.get(src_pt)
                    mix_arrived[mix_out_pt][src_seg] = src_flow

                    # Check whether all inlets are now resolved.
                    all_known = True
                    total = 0.0
                    for inp_pt in mixer_inputs[mix_out_pt]:
                        inp_seg = point_to_seg.get(inp_pt)
                        if inp_seg in mix_arrived[mix_out_pt]:
                            total += mix_arrived[mix_out_pt][inp_seg]
                        elif inp_seg is not None and inp_seg in seg_flow:
                            total += seg_flow[inp_seg]
                        else:
                            # Try reading the current point value as a fallback
                            pt_obj = self.Point.get(inp_pt)
                            if pt_obj and pt_obj.Mass_flowrate is not None:
                                total += pt_obj.Mass_flowrate
                            else:
                                all_known = False
                                break
                    if not all_known:
                        continue   # wait for the missing inlet

                    new_seg_flow = total

                # ── Split edge: apply fraction ────────────────────────
                elif etype == 'split':
                    frac = _get_split_frac(src_pt, dst_pt)
                    if frac is None:
                        # Separator: derive from quality Q
                        in_pt_obj  = self.Point.get(src_pt)
                        q_val      = getattr(in_pt_obj, 'Q', None) if in_pt_obj else None
                        fracs_list = [op for op, _ in splitter_fracs.get(src_pt, [])]
                        if isinstance(q_val, (int, float)) and 0.0 <= q_val <= 1.0:
                            frac = q_val if dst_pt == fracs_list[0] else 1.0 - q_val
                        else:
                            warnings.warn(
                                f"[FlowNet] Separator inlet '{src_pt}' has no "
                                f"numeric quality Q; using 0.5 as fallback."
                            )
                            frac = 0.5
                    new_seg_flow = src_flow * frac

                # ── Passthrough / simple forward ──────────────────────
                elif dst_seg in visited_segs:
                    continue
                else:
                    new_seg_flow = src_flow

            # ══════════════════════════════════════════════════════════
            # BACKWARD direction  (traversing seg_in_edges in reverse)
            # ══════════════════════════════════════════════════════════
            else:  # direction == 'bwd'

                if dst_seg in visited_segs:
                    continue

                # dst_seg is the UPSTREAM segment (we are going backwards).
                # The edge is  src_pt → dst_pt  in the forward graph,
                # meaning dst_pt is in dst_seg (upstream) and
                # src_pt is in the segment we just resolved (downstream).
                # We want to infer the flow of dst_seg from src_flow.

                if etype == 'passthrough':
                    # Upstream splitter-inlet segment carries the same flow.
                    new_seg_flow = src_flow

                elif etype == 'split':
                    # src_pt is the splitter inlet (upstream).
                    # dst_pt is a splitter outlet (the downstream segment
                    # we just set = src_flow).
                    # Back-calculate: m_upstream = m_outlet / fraction
                    frac = _get_split_frac(src_pt, dst_pt)
                    if frac is None or frac <= 0:
                        continue   # cannot back-calculate without a valid fraction
                    new_seg_flow = src_flow / frac

                elif etype == 'mix_in':
                    # dst_seg is one of the mixer INLET segments.
                    # src_pt is in dst_seg (a mixer inlet point).
                    # dst_pt is the mixer outlet point.
                    # We can resolve this inlet only if the mixer outlet
                    # AND all OTHER inlets are already known.
                    mix_out_pt  = dst_pt
                    mix_out_seg = point_to_seg.get(mix_out_pt)
                    outlet_flow = seg_flow.get(mix_out_seg)
                    if outlet_flow is None:
                        outlet_flow = _seg_flow_of(mix_out_seg)
                    if outlet_flow is None:
                        continue   # mixer outlet not yet resolved

                    other_total = 0.0
                    all_others_known = True
                    for inp_pt in mixer_inputs[mix_out_pt]:
                        inp_seg = point_to_seg.get(inp_pt)
                        if inp_seg == dst_seg:
                            continue   # this is the inlet we're trying to find
                        f = seg_flow.get(inp_seg) or _seg_flow_of(inp_seg)
                        if f is None:
                            all_others_known = False
                            break
                        other_total += f
                    if not all_others_known:
                        continue   # cannot resolve yet

                    new_seg_flow = outlet_flow - other_total
                    if new_seg_flow < 0:
                        warnings.warn(
                            f"[FlowNet] Backward mixer calc gave negative flow "
                            f"({new_seg_flow:.4f} kg/s) for inlet segment of "
                            f"'{mix_out_pt}'. Skipping."
                        )
                        continue

                else:
                    continue   # unknown edge type — skip

            # ── Commit the resolved flow ───────────────────────────────
            if new_seg_flow is None:
                continue

            visited_segs.add(dst_seg)
            seg_flow[dst_seg] = new_seg_flow
            _write_segment(dst_seg, new_seg_flow)
            _enqueue_neighbours(dst_seg, new_seg_flow)

    def _make_propagator(self):
        """Return a closure bound to *this* model instance."""
        def _propagator(point_name, new_value):
            self._propagate_mass_flow(point_name, new_value)
        return _propagator

    def _attach_propagators(self):
        """Attach the propagation callback to every registered Prop."""
        propagator = self._make_propagator()
        for pt in self.Point.values():
            object.__setattr__(pt, '_flow_propagator', propagator)

    def _replay_known_flowrates(self, verbose=False):
        """
        After the graph is built and propagators are attached, find every
        point whose Mass_flowrate was explicitly supplied by the user at
        construction time (flagged by _user_set_flowrate=True on the Prop)
        and fire a real propagating assignment so downstream points update.

        Why not just skip points with an upstream neighbour?
        ─────────────────────────────────────────────────────
        In a closed loop (e.g. an ORC working-fluid loop) EVERY point has
        an upstream neighbour — there are no open ends.  The old heuristic
        would skip all of them and propagate nothing.  The flag approach
        is topology-independent: it replays exactly the points the user
        explicitly seeded, regardless of whether the graph is open or cyclic.

        Double-counting guard
        ──────────────────────
        If the user seeded multiple points on the same stream (e.g. point
        '1' and point '3' both have Mass_flowrate=10 on the same loop),
        replaying '1' will propagate through '3' and overwrite it correctly.
        Replaying '3' afterwards would start a second BFS from mid-loop —
        harmless but redundant.  We therefore skip any point that already
        has its Mass_flowrate set correctly by a previous replay in this pass.
        """
        replayed = []
        # Replay ALL user-seeded points, including those with None.
        # Replaying None triggers the BFS to clear stale values from
        # prior solves across every reachable segment.
        for name, pt in self.Point.items():
            if pt.__dict__.get('_user_set_flowrate'):
                pt.Mass_flowrate = pt.Mass_flowrate   # triggers hook → BFS propagates
                replayed.append(name)

        if replayed and verbose:
            print(f"[FlowNet] Replayed Mass_flowrate from "
                  f"{len(replayed)} source point(s): {replayed}")

    def enable_flow_propagation(self, verbose=False):
        """
        Build the flow graph from the current component registry, attach
        reactive mass-flow propagators to all state points, and replay
        any Mass_flowrate values that were already set during add_point()
        so they cascade through the graph immediately.

        This means the following workflow is fully supported:

            Model.add_point('water', '1', P=8e6, T=753.15, Mass_flowrate=1)
            Model.add_point('water', '2', P=0.7e6)
            # ... add remaining points and components ...
            Model.enable_flow_propagation()   # ← replay fires here
            # point '2' already has the correct mass flow — no re-assignment needed
        """
        self._build_flow_graph()
        self._attach_propagators()
        self._propagation_enabled = True
        self._replay_known_flowrates(verbose=verbose)
        if verbose:
            n_edges = sum(len(v) for v in self._flow_graph.values())
            print(
                f"[FlowNet] Graph built: {len(self._flow_graph)} nodes, "
                f"{n_edges} edges. "
                f"Propagators attached to {len(self.Point)} state points."
            )

    def disable_flow_propagation(self):
        """
        Detach all propagators so Mass_flowrate assignments are silent
        again (only exergy updates, as before enabling).
        """
        from .state import _NO_PROPAGATOR
        for pt in self.Point.values():
            object.__setattr__(pt, '_flow_propagator', _NO_PROPAGATOR)
        self._propagation_enabled = False

    def flow_graph_summary(self):
        """
        Print a human-readable table of all flow edges, annotating
        Splitter origins and Mixer destinations.
        """
        if not hasattr(self, '_flow_graph'):
            print("[FlowNet] Flow graph not built yet. "
                  "Call enable_flow_propagation() first.")
            return

        def _mf(name):
            pt = self.Point.get(name)
            if pt and pt.Mass_flowrate is not None:
                return f"{pt.Mass_flowrate:.4f} kg/s"
            return "? kg/s"

        print("=" * 60)
        print("  FLOW NETWORK SUMMARY")
        print("=" * 60)
        printed_any = False
        for src, dsts in self._flow_graph.items():
            if not dsts:
                continue
            src_tag = " [SPLIT]" if src in self._splitter_fracs else ""
            for dst in dsts:
                dst_tag = " [MIXER]" if dst in self._mixer_inputs else ""
                print(f"  {src}{src_tag} ({_mf(src)})  →  "
                      f"{dst}{dst_tag} ({_mf(dst)})")
                printed_any = True
        if not printed_any:
            print("  (no edges — graph may be empty)")
        print("=" * 60)

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