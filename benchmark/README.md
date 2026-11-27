# Benchmark and supplement scripts

The scripts and raw results behind Section 2.4 of the article and its supplement:
the comparison with TESPy and the further checks of the solver.

Run each script from its own folder with ThermoSim installed (`sweep.py` takes the
number of unknowns: `python sweep.py 2`, `python sweep.py 3`). The versions used
are in `../requirements-tested.txt` (TESPy is needed only here). Several scripts
write a JSON file into the current folder, replacing the committed result.
Timings depend on the machine.

| script | what it does |
|---|---|
| p1_thermosim.py, p1_tespy.py | geothermal evaporator, pinch specified (S9, Table S6) |
| p234_thermosim.py, p234_tespy.py | crossing, unsolvable and hot-side-unknown specifications (S9, Table S7) |
| near_miss.py, case8_9.py, cases12.py | approach below the requested value; cases 8, 9 and 1, 2 |
| sweep.py, sweep3_condenser.py, sweep3_consistent.py | every two- and three-unknown specification (`sweep_*.json`, `*.log`) |
| timing.py, timing_vs_divisions.py | time per solve, and against the number of divisions (S9, Table S7d) |
| thr.py | width of the exclusion window (S15) |
| hrsg.py, feas.py | the combined-cycle designs (S8, S13) |
| orc.py, exergy_grid.py | exergy balance and the pressure and approach grid (S11, S12) |
| disc.py, dp_rb.py, rb2.py | discretisation, double-pinch and rollback checks |
| flow.py | dispatch flowchart (`flow.png`, Figure S1) |

Further scripts are described in [round2/README.md](round2/README.md).

## Notes on reproducing

- The committed JSON and log files are the original results. The tables of the
  article for all 35 specifications come from `round2/m3.py` and
  `round2/tespy_sweep.py`.
- `sweep.py` is an earlier sweep on a state whose cold outlet is saturated
  vapour, which is not consistent with a hot-end anchor. For some three-unknown
  cases it therefore returns a different design, and `{H_hi, H_ho, mh}` can
  return a physically meaningless flow rate. `sweep3_consistent.py` repeats
  those cases on a consistent state.
- TESPy outcomes can differ between platforms. On macOS (arm64), TESPy 0.11.3
  returns a different state for the three-unknown case `{H_hi, H_ho, H_co}`
  where the committed result reports no solution, giving 5 and 10 instead of
  4 and 11. The ThermoSim outcomes are the same.
