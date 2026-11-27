# Further scripts for the supplement

Run from this directory with the package on the path, e.g.
`PYTHONPATH=/path/to/ThermoSim python script.py` (tested: ThermoSim 3.5.0, CoolProp 8.0.0,
SciPy 1.18.1, NumPy 2.5.3, TESPy 0.11.3, Python 3.13).

| script | supplement section |
|---|---|
| m3.py, tespy_sweep.py | S7b, S7c: all 35 specifications, ThermoSim and TESPy |
| ../timing_vs_divisions.py, ../timing.py | S7d: time per solve |
| hrsg3.py | S4b: three-unknown heat-recovery steam generator |
| c345.py | S14: family of designs, Jacobian, same-side example (S16) |
| conv.py, ../thr.py | S15: discretisation, exclusion window |
| chena_s.py, chena_sens.py, chena_prop.py, chena_pc.py | S18: Chena sensitivity, equation of state, round trip |
| elas.py | S19: elasticities |
| rbk.py | state rollback checks |

chena_s.py imports ThermoSim only; the *.json files hold the results quoted in the
supplement. ../thr.py and tespy_sweep.py write their JSON next to themselves.
