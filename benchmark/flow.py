import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams['font.family']='DejaVu Sans'
fig,ax=plt.subplots(figsize=(6.6,9.2)); ax.set_xlim(0,100); ax.set_ylim(0,150); ax.axis('off')
B='#1f4e79'; R='#9c2a1f'; G='#2e7d32'; FILL='#eef3f8'
def box(x,y,w,h,t,fc=FILL,ec=B,fs=7.4,bold=False):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.2,rounding_size=1.2',fc=fc,ec=ec,lw=1.1))
    ax.text(x,y,t,ha='center',va='center',fontsize=fs,fontweight='bold' if bold else 'normal',linespacing=1.3)
def arr(x1,y1,x2,y2,c=B,ls='-'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=8,lw=1.0,color=c,ls=ls,shrinkA=0,shrinkB=0))
def lab(x,y,t,fs=7.4,c=B,b=True): ax.text(x,y,t,ha='center',va='center',fontsize=fs,color=c,fontweight='bold' if b else 'normal')
box(50,146,68,6,'Count the unknowns n among H$_{hi}$, H$_{ho}$, H$_{ci}$, H$_{co}$, $\\dot{m}_h$, $\\dot{m}_c$',bold=True)
box(50,136,68,5,'Capture the four state points and both mass flow rates')
arr(50,143,50,138.7)
xs=[13,38,62,87]; ttl=['n = 0','n = 1','n = 2','n = 3']
txt=['Verify that the\ntwo sides agree\non the duty;\nbuild the profile','The energy balance\ngives the unknown\ndirectly','Solve the pair\nfor the pinch\n(panel a)','Anchor one unknown,\nthen solve the\nremaining pair\n(panel b, then a)']
ax.plot([13,87],[127.5,127.5],color=B,lw=1); arr(50,133.3,50,127.5) if False else ax.plot([50,50],[133.4,127.5],color=B,lw=1)
for x,t,tx in zip(xs,ttl,txt):
    arr(x,127.5,x,123.2); lab(x,125.6,t) if False else None
    ax.text(x+0.0,129.6,t,ha='center',va='center',fontsize=7.4,color=B,fontweight='bold')
    box(x,116,22,13,tx,fs=7.0)
box(50,98,80,14,'Validate the solved exchanger:\n(1) the two sides agree on the duty;  (2) terminal approach ≥ ΔT$_{pp}$;\n(3) whole-profile minimum approach ≥ ΔT$_{pp}$.\nA violation is reported with the computed minimum: warning,\nor error if strict_pinch is set',fs=7.0)
for x in xs: ax.plot([x,x],[109.5,105.5],color=B,lw=1)
ax.plot([13,87],[105.5,105.5],color=B,lw=1); arr(50,105.5,50,104.7)
box(50,84,62,6,'Solved: state points written, profile available',fc='#e8f3ea',ec=G,bold=True,fs=7.2)
arr(50,91.5,50,87.2)
# panels
ax.plot([4,96],[75.5,75.5],color='#999999',lw=0.8,ls=':')
lab(26,72.5,'(a) Pair solve (n = 2)',fs=7.6)
lab(74,72.5,'(b) Anchor (n = 3)',fs=7.6)
ya=[64,52,40,28]
ta=['Choose the free variable x\nfrom Table S1','Bracket x from the physical floor\n(no state is colder than the cold inlet)','Brent: solve  f(x) = min(T$_h$ − T$_c$) − ΔT$_{pp}$ = 0\nNo sign change: narrower bracket,\nthen differential evolution','Root found with |f| small?\nNo: refuse and name the cause']
for y,t in zip(ya,ta): box(26,y,44,9.5,t,fs=6.9)
for i in range(3): arr(26,ya[i]-4.8,26,ya[i+1]+4.8)
tb=['(i) Is the anchor admissible\nfor this exchanger type?','(ii) Does it reduce to a solvable\n2-unknown case?','(iii) Do the temperatures stay\nordered (T$_{hi}$ > T$_{ho}$, T$_{ci}$ < T$_{co}$)?','Apply the anchor and exclude its end\nfrom the search; solve the remaining\npair as in panel (a)']
for y,t in zip(ya,tb): box(74,y,44,9.5,t,fs=6.9)
for i in range(3): arr(74,ya[i]-4.8,74,ya[i+1]+4.8)
lab(74,19.2,'A check fails: try the next anchor; none left: refuse,\nmodel unchanged',fs=6.5,c=R,b=False)
box(50,9,92,12,'Any exception, whether from the anchor checks, the root-finder or the validation:\nthe captured state points and mass flow rates are restored before the error is raised.\nThe error names the cause (n > 3, unsolvable combination,\nno root found, or approach violated under strict_pinch).',fc='#fbeeee',ec=R,fs=6.9)
plt.savefig('flow.png',dpi=240,bbox_inches='tight',facecolor='white'); print('ok')
