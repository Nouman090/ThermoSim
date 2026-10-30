# `nh3h2o.py` — a CoolProp-style NH₃–H₂O property library

Give it the overall ammonia mass fraction plus **any valid pair** of properties and
it returns the complete state.

```python
import nh3h2o as aw

st = aw.state(z=0.4, T=350, P=10)      # full state object
st = aw.state(z=0.4, P=10, Q=0)        # bubble point
st = aw.state(z=0.4, P=10, H=500)      # PH flash
st.T, st.Q, st.x, st.y, st.h, st.s, st.v, st.u, st.rho, st.phase

aw.PropsSI('T', 'P', 10, 'Q', 0, 0.4)  # single value, CoolProp style
aw.bubble_T(10, 0.4);  aw.dew_T(10, 0.4)
aw.bubble_T(10, 0.4, vle='PK')         # Patek–Klomfar instead
```

**Accepted inputs** — `T`, `P`, `H`, `S`, `U`, `V`, `D`, `Q` (case-insensitive; `D`
is density, `V` specific volume). Exactly two, plus `z`. One of them must be `T` or `P`.

**Units** — T [K], P [bar], h [kJ/kg], s [kJ/kg·K], v [m³/kg], u [kJ/kg]. All
compositions are ammonia **mass** fractions. Reference: h = s = 0 for the saturated
liquid of each pure component at 273.15 K.

## Model

| Layer | Source |
|---|---|
| Single-phase properties | Ziegler & Trepp (1984), Gibbs free energy |
| Vapour–liquid equilibrium | solved from the **same** Gibbs function (equal chemical potentials) |
| Optional VLE | Pátek & Klomfar (1995), via `vle='PK'` |

Solving equilibrium from the same fundamental relation that supplies the enthalpies is
what makes energy balances close exactly. Pairing Ziegler–Trepp enthalpies with the
Pátek–Klomfar correlations — the common shortcut — leaves a 1–5 K bubble-point
mismatch below 20 bar and up to 11 K at 50 bar, which shows up as a residual in every
component balance. `vle='PK'` is provided for comparison, not as the default.

## How the flash works

For a Gibbs-excess model the equilibrium conditions rearrange to

```
K_A = y/x       = exp{ [g_A^L − g_A^V + G^E + (1−x)·dG^E/dx] / RT }
K_W = (1−y)/(1−x) = exp{ [g_W^L − g_W^V + G^E − x·dG^E/dx] / RT }
```

Both depend on **(T, P, x) only — never on y**. So the bubble surface is the single
scalar equation `K_A·x + K_W·(1−x) = 1`, and every flash reduces to nested 1-D root
finds with guaranteed brackets (the two-phase region always lies between the two
pure-component saturation temperatures) instead of a fragile multivariate solve.
`dG^E/dx` is differentiated analytically at import, so the chemical potentials are exact.

## Validation — 22/22 checks pass

Reference data: CoolProp 8.0.0 (IAPWS-95 for water, Tillner-Roth for ammonia).
Full log in `validation_output.txt`; run `python validate_lib.py` to reproduce.

| Test | Result |
|---|---|
| Enthalpy vs reference, pure limits, both phases, P ≤ 40 bar | AAD 1.38 kJ/kg, max 5.53 kJ/kg |
| Entropy, max deviation | 0.015 kJ/kg·K |
| Specific volume, max deviation | 0.94% |
| Bubble & dew temperature at pure limits | max error 0.27 K |
| Latent heat, max deviation | 1.22% |
| Round trip, all 9 input pairs (80 base states) | inputs reproduced to 2×10⁻⁸% |
| Two-phase mass, species, energy, entropy, volume balances | ≤ 2×10⁻¹⁶ |
| K-value closure (VLE consistency) | 9×10⁻¹³ |
| Enthalpy continuity across bubble and dew lines | 5×10⁻⁴ kJ/kg |
| Absorption-cycle energy balance (`demo_cycle.py`) | 6.6×10⁻¹² % of Q_gen |

Performance: **0.18 ms** per TP flash, 16 ms per PH flash.

## Non-uniqueness — read this before using `T,H` or `P,V`

Two input pairs are genuinely multi-valued, and the library **warns** rather than
silently picking a branch:

- **`(T, H)` in the compressed liquid.** Liquid enthalpy is nearly independent of
  pressure, so many pressures satisfy the same enthalpy. 35 of 80 test states were
  ambiguous. If the residual is flat across the whole range the library raises
  `AmmoniaWaterError` telling you the inputs do not determine the state.
- **`(P, V)` in the compressed liquid.** Liquid v(T) passes through a minimum (the
  density extremum), so two temperatures give the same volume. 11 of 80 were ambiguous.

Every root returned does satisfy the inputs — they are alternative valid states, not
errors. For cycle work use `(P,T)`, `(P,H)`, `(P,S)` or `(P,Q)`, all of which were
single-valued at every test point.

## Validity

**230–500 K, 0.2–50 bar** (Ziegler & Trepp's validated envelope). Outside it the
library warns. Accuracy degrades sharply above 50 bar — the vapour Gibbs function
contains C₄·P_r³/T_r¹¹ terms, so error grows with the cube of pressure and reaches
−36% at 100 bar. For pressures above 50 bar use Ibrahim & Klein (1993, to 110 bar)
or Tillner-Roth & Friend (1998).

Pátek & Klomfar (`vle='PK'`) is itself only fitted to 2 MPa (20 bar) and is not
thermodynamically consistent — its bubble and dew surfaces can intersect.

## Files

| File | Purpose |
|---|---|
| `nh3h2o.py` | the library |
| `validate_lib.py` | 22-check validation suite |
| `demo_cycle.py` | worked absorption chiller with balance checks |
| `validation_output.txt` | full log of both |
| `fig8_library_parity.png` | parity of h, s, v vs reference |
| `fig9_ph_ts_diagrams.png` | P–h and T–s diagrams built through the flash API |

## References

- Ziegler, B. & Trepp, Ch. (1984). *Int. J. Refrigeration* **7**(2), 101–106.
- Pátek, J. & Klomfar, J. (1995). *Int. J. Refrigeration* **18**(4), 228–234.
- Ibrahim, O.M. & Klein, S.A. (1993). *ASHRAE Trans.* **99**(1), 1495–1502.
- Tillner-Roth, R. & Friend, D.G. (1998). *J. Phys. Chem. Ref. Data* **27**(1), 63–96.
