"""
Figure -- Chena Hot Springs geothermal ORC.

(a) Deviation of ThermoSim and of the published IPSEpro simulation from the
    measured plant data.  Two ThermoSim runs are shown: at the plant's own
    inputs (16.00 bar, brine exit 54.44 C) and at the inputs IPSEpro used
    (16.95 bar, brine exit 54.94 C), so that the two codes are compared at
    common inputs.  The cooling-water outlet is compared on the kelvin scale
    (percentages of a Celsius reading depend on the temperature origin); the
    absolute difference in K is printed under the bar.
(b) Evaporator profiles at the three operating points that matter.

Regenerates fig_chena.png.   Run with:  MPLBACKEND=Agg python fig_chena.py
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, warnings, pandas as pd
warnings.simplefilter('ignore')
import validate_chena as V
from ThermoSim.plotting import CyclePlotter

def results(**kw):
    M = V.build(**kw)
    return M, dict(m=M.Point['1'].Mass_flowrate, cw=M.Point['c2'].T,
                   Qe=M.Component['Evaporator'].Q/1e3, Qc=M.Component['Condenser'].Q/1e3)

M_plant, r_plant = results()                                       # 16.00 bar, 54.44 C
M_fit,   r_fit   = results(P_evap=16.95e5)                         # 16.95 bar, 54.44 C
M_ips,   r_ips   = results(P_evap=16.95e5, T_geo_out_C=54.94)      # IPSEpro inputs

labels = ['working fluid\nflow rate', 'cooling water\noutlet (K scale)', 'evaporator\nduty', 'condenser\nduty']
plant = dict(m=12.17, cw=10.00 + 273.15, Qe=2580.0, Qc=2360.0)
ipse  = dict(m=11.99, cw=9.91 + 273.15,  Qe=2570.38, Qc=2327.10)
keys = ['m', 'cw', 'Qe', 'Qc']
dev = lambda r: [100*(r[k]-plant[k])/plant[k] for k in keys]
dK  = lambda r: r['cw'] - plant['cw']
series = [('ThermoSim, plant inputs (16.00 bar)', dev({**r_plant, 'cw': r_plant['cw']}), dK(r_plant), '#1f6fb4'),
          ('ThermoSim, IPSEpro inputs (16.95 bar)', dev(r_ips), dK(r_ips), '#2e7d32'),
          ('IPSEpro [24]', dev(ipse), dK(ipse), '#c0392b')]

fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1), gridspec_kw={'width_ratios': [1.35, 1]})
ax = axes[0]
x = np.arange(len(labels)); w = 0.26
for j, (name, e, dk, col) in enumerate(series):
    xs = x + (j-1)*w
    ax.bar(xs, e, w, label=name, color=col)
    for i, (xi, v) in enumerate(zip(xs, e)):
        txt = f'{v:+.2f}'
        ax.text(xi, v + (0.10 if v >= 0 else -0.10), txt, ha='center',
                va='bottom' if v >= 0 else 'top', fontsize=6.6, color=col, fontweight='bold')
    ax.text(x[1] + (j-1)*w, 0.55 + 0.0*j, f'{dk:+.2f} K', ha='center', va='bottom',
            fontsize=6.4, color=col, fontweight='bold', rotation=90)
ax.axhline(0, color='#333', lw=1.0)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.2)
ax.set_ylabel('signed deviation from plant data [%]\n(model - plant)/plant')
ax.set_title('(a)  Chena geothermal ORC', fontsize=10)
ax.legend(frameon=False, fontsize=7.6, loc='upper left'); ax.grid(axis='y', alpha=0.25); ax.set_axisbelow(True)
ax.set_ylim(-2.6, 3.4)

ax = axes[1]
pts = [('16.00 bar, brine 54.44 °C', M_plant, '-', '#1f6fb4'),
       ('16.95 bar, brine 54.44 °C', M_fit, '--', '#7b3fa0'),
       ('16.95 bar, brine 54.94 °C (IPSEpro)', M_ips, ':', '#2e7d32')]
first = True; mins = []
for k, (tag, Mx, ls, col) in enumerate(pts):
    p = CyclePlotter(Mx)
    p.plot_hex_profile('Evaporator', div_N=200, save_csv=True, csv_path=f'/tmp/ch{k}.csv')
    d = pd.read_csv(f'/tmp/ch{k}.csv'); Q = d['Heat_Transferred_kW']/1e3
    if first:
        ax.plot(Q, d['Hot_Temperature_C'], '-', color='#444', lw=2.0, label='geothermal brine'); first = False
    mn_ = d['Delta_T_K'].min()
    ax.plot(Q, d['Cold_Temperature_C'], ls, color=col, lw=1.8, label=f'R134a, {tag}: min ΔT = {mn_:.2f} K')
    mins.append((tag, mn_, col))
ax.set_xlabel('cumulative duty [MW]'); ax.set_ylabel('temperature [°C]')
ax.set_title('(b)  evaporator profile, three operating points', fontsize=10)
ax.grid(alpha=0.25); ax.set_ylim(5, 85)
ax.legend(frameon=False, fontsize=6.8, loc='lower right')
fig.tight_layout()
fig.savefig('fig_chena.png', dpi=300, bbox_inches='tight')
for tag, mn, _ in mins: print(tag, '->', round(mn, 3), 'K')
for n, e, dk, _ in series: print(n, [round(v, 3) for v in e], 'cw dK', round(dk, 3))
