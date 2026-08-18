import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from validate_moran import reheat

etas = np.linspace(0.85, 1.00, 31)
res = []
for e in etas:
    M, _ = reheat(eta_t=float(e), eta_p=1.0)
    M.ModelSummary(verbose=False)
    res.append(M.Efficiency/100.0)

print("eta_t   cycle eta")
for e, r in zip(etas[::6], res[::6]):
    print(f" {e:.2f}     {r:.4f}")
print(f"\nendpoint check:")
print(f"  eta_t=0.85 -> {res[0]:.4f}   (Moran Ex 8.4 = 0.351)")
print(f"  eta_t=1.00 -> {res[-1]:.4f}   (Moran Ex 8.3 = 0.403)")

fig, ax = plt.subplots(figsize=(6.2, 4.8))
ax.plot(etas, res, '-', color='#2e7d32', lw=2.2, label='ThermoSim')
ax.plot([0.85, 1.00], [0.351, 0.403], 'o', ms=8, mfc='none', mec='#c0392b',
        mew=1.8, label='Moran Ex. 8.4 / 8.3')
ax.set_xlabel('Isentropic turbine efficiency')
ax.set_ylabel('Cycle thermal efficiency')
ax.set_xlim(0.85, 1.00); ax.set_ylim(0.32, 0.42)
ax.set_yticks(np.arange(0.32, 0.43, 0.02))
ax.grid(alpha=0.25); ax.legend(frameon=False, loc='upper left')
fig.tight_layout()
fig.savefig('/tmp/moran_fig.png', dpi=110)
print("figure saved")
