import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, warnings, pandas as pd
warnings.simplefilter('ignore')
import validate_chena as V
from ThermoSim.plotting import CyclePlotter

M = V.build()
ts = dict(m=M.Point['1'].Mass_flowrate,
          cw=M.Point['c2'].T-273.15,
          Qe=M.Component['Evaporator'].Q/1e3,
          Qc=M.Component['Condenser'].Q/1e3)
labels = ['working fluid\nflow rate', 'cooling water\noutlet', 'evaporator\nduty', 'condenser\nduty']
plant  = [12.17, 10.00, 2580.0, 2360.0]
ips    = [11.99,  9.91, 2570.38, 2327.10]
mine   = [ts['m'], ts['cw'], ts['Qe'], ts['Qc']]
e_ts = [100*(a-b)/b for a,b in zip(mine, plant)]
e_ip = [100*(a-b)/b for a,b in zip(ips,  plant)]

fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9),
                         gridspec_kw={'width_ratios':[1.15,1]})

# (a) deviation from measured plant data
ax = axes[0]
x = np.arange(len(labels)); w = 0.36
ax.bar(x-w/2, e_ts, w, label='ThermoSim', color='#1f6fb4')
ax.bar(x+w/2, e_ip, w, label='IPSEpro (Ref. [11])', color='#c0392b')
ax.axhline(0, color='#333', lw=1.0)
for xi, v in zip(x-w/2, e_ts):
    ax.text(xi, v + (0.12 if v>=0 else -0.28), f'{v:+.2f}', ha='center',
            fontsize=8, color='#1f6fb4', fontweight='bold')
for xi, v in zip(x+w/2, e_ip):
    ax.text(xi, v + (0.12 if v>=0 else -0.28), f'{v:+.2f}', ha='center',
            fontsize=8, color='#c0392b', fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylabel('deviation from measured plant data [%]')
ax.set_title('(a)  Chena geothermal ORC', fontsize=10)
ax.legend(frameon=False, fontsize=8.5, loc='upper left')
ax.grid(axis='y', alpha=0.25); ax.set_axisbelow(True)
ax.set_ylim(-2, 3)

# (b) evaporator pinch at the two candidate turbine inlet pressures
ax = axes[1]
res = []
for P, tag, col in ((16.00e5,'plant  16.00 bar','#1f6fb4'),
                    (16.95e5,'fitted 16.95 bar (Ref. [11]])','#c0392b')):
    V.P_EVAP = P
    Mx = V.build()
    p = CyclePlotter(Mx)
    p.plot_hex_profile('Evaporator', div_N=200, save_csv=True,
                       csv_path=f'/tmp/ch{int(P)}.csv')
    d = pd.read_csv(f'/tmp/ch{int(P)}.csv')
    Q = d['Heat_Transferred_kW']/1e3
    if not res:                       # brine stream is the same in both cases
        ax.plot(Q, d['Hot_Temperature_C'], '-', color='#444', lw=2.0,
                label='geothermal brine')
        ax.plot(Q, d['Cold_Temperature_C'], '-', color=col, lw=2.0,
                    label=f'R134a, {tag.split()[1]} bar')
    else:
        ax.plot(Q, d['Cold_Temperature_C'], '--', color=col, lw=1.5,
                label=f'R134a, {tag.split()[1]} bar')
    res.append((tag, d['Delta_T_K'].min(), col))
for i,(tag, mn, col) in enumerate(res):
    ax.text(0.035, 0.93-0.085*i, f'{tag}:  min $\\Delta T$ = {mn:.2f} K',
            transform=ax.transAxes, fontsize=8.6, color=col, fontweight='bold')
ax.set_xlabel('cumulative duty [MW]'); ax.set_ylabel('temperature [°C]')
ax.set_title('(b)  evaporator profile, two reported pressures', fontsize=10)
ax.legend(frameon=False, fontsize=8.2, loc='lower right'); ax.grid(alpha=0.25)
ax.set_ylim(5, 82)

fig.tight_layout()
fig.savefig('fig_chena.png', dpi=300, bbox_inches='tight')
for tag, mn, _ in res: print(tag, '->', round(mn,3), 'K')
