"""
Figure -- the anchoring rule.
(a) The four relations by which a known temperature at one terminal fixes the
    temperature facing it across the approach.  The cold-end relations are
    placed at the cold end and the hot-end relations at the hot end.
(b) Terminals at which an anchor is admissible, by exchanger type.
Regenerates fig_anchor_rules.png.   Run with:  MPLBACKEND=Agg python fig_anchor_rules.py
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
RED, BLUE, BRN, SHADE = '#b5381e', '#1f5fa0', '#8a5a00', '#f7dfb4'
fig = plt.figure(figsize=(12.4, 5.2))
gs = GridSpec(3, 2, width_ratios=[1.55, 1], hspace=0.55, wspace=0.18, figure=fig)
ax = fig.add_subplot(gs[:, 0])
x0, x1 = 0.17, 0.83
Th = (0.30, 0.90); Tc = (0.15, 0.78)          # (cold end, hot end) temperatures, schematic
ax.plot([x0, x1], Th, color=RED, lw=3.2); ax.plot([x0, x1], Tc, color=BLUE, lw=3.2, ls=(0, (5, 2.2)))
for x, y, c in ((x0, Th[0], RED), (x1, Th[1], RED), (x0, Tc[0], BLUE), (x1, Tc[1], BLUE)): ax.plot(x, y, 'o', color=c, ms=7)
for x in (x0, x1):
    lo, hi = (Tc[0], Th[0]) if x == x0 else (Tc[1], Th[1])
    ax.annotate('', xy=(x, lo), xytext=(x, hi), arrowprops=dict(arrowstyle='<->', color='#555', lw=1.1))
ax.text(x0-0.03, Th[0]+0.02, '$T_{ho}$', ha='right', fontsize=14); ax.text(x0-0.03, Tc[0]-0.03, '$T_{ci}$', ha='right', fontsize=14)
ax.text(x1+0.03, Th[1]+0.01, '$T_{hi}$', ha='left', fontsize=14); ax.text(x1+0.03, Tc[1]-0.03, '$T_{co}$', ha='left', fontsize=14)
ax.text(x0+0.20, 0.235, 'at the cold end', color=BRN, fontsize=13)
ax.text(x0+0.20, 0.165, '$T_{ci} = T_{ho} - \\Delta T_{pp}$', fontsize=14)
ax.text(x0+0.20, 0.095, '$T_{ho} = T_{ci} + \\Delta T_{pp}$', fontsize=14)
ax.text(x1-0.02, 0.43, 'at the hot end', color=BRN, fontsize=13, ha='right')
ax.text(x1-0.02, 0.365, '$T_{co} = T_{hi} - \\Delta T_{pp}$', fontsize=14, ha='right')
ax.text(x1-0.02, 0.30, '$T_{hi} = T_{co} + \\Delta T_{pp}$', fontsize=14, ha='right')
ax.text(x0, 0.07, 'cold end', ha='center', color='#555', fontsize=12); ax.text(x1, 0.07, 'hot end', ha='center', color='#555', fontsize=12)
ax.set_xlim(0, 1); ax.set_ylim(0.04, 1.0)
ax.set_xlabel('cumulative heat-transfer rate (duty) $\\dot Q$ →', color='#555', fontsize=12); ax.set_ylabel('temperature', color='#555', fontsize=12)
ax.set_xticks([]); ax.set_yticks([])
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('(a)  the four anchor relations', loc='left', fontsize=13)
def mini(pos, title, hot, cold, ends):
    a = fig.add_subplot(gs[pos, 1]); xs = [0, 0.1, 0.3, 0.7, 0.9, 1.0]
    a.plot(xs[:len(hot)] if len(hot) != 6 else xs, hot, color=RED, lw=2.4); a.plot(xs, cold, color=BLUE, lw=2.4, ls=(0, (4, 2)))
    for e in ends: a.axvspan(*( (0, 0.1) if e == 'cold' else (0.9, 1.0) ), color=SHADE)
    a.set_title(title, fontsize=11); a.set_xticks([]); a.set_yticks([])
    for s in ('top', 'right'): a.spines[s].set_visible(False)
    return a
ax1 = mini(0, 'counter-flow: either end', [0.40, 0.45, 0.55, 0.75, 0.90, 0.95], [0.20, 0.28, 0.38, 0.62, 0.80, 0.88], ['cold', 'hot'])
ax2 = mini(1, 'evaporator: hot-inlet end', [0.40, 0.45, 0.55, 0.75, 0.90, 0.95], [0.20, 0.45, 0.45, 0.45, 0.45, 0.80], ['hot'])
ax3 = mini(2, 'condenser: hot-outlet end', [0.45, 0.60, 0.60, 0.60, 0.60, 0.85], [0.20, 0.30, 0.42, 0.58, 0.70, 0.78], ['cold'])
ax1.set_title('(b)  where an anchor is admissible\ncounter-flow: either end', fontsize=11)
fig.savefig('fig_anchor_rules.png', dpi=250, bbox_inches='tight')
