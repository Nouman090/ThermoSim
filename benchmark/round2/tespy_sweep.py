"""All 15 two-unknown and 20 three-unknown specifications of the benchmark evaporator in TESPy 0.11.3
(MovingBoundaryHeatExchanger, td_pinch = 8 K; ttd_u = 8 K added for three unknowns)."""
import itertools, time, json, logging, io, warnings, sys
warnings.filterwarnings('ignore')
logging.disable(logging.CRITICAL)
import numpy as np
from CoolProp.CoolProp import PropsSI as P
from tespy.networks import Network
from tespy.components import Source, Sink, MovingBoundaryHeatExchanger
from tespy.connections import Connection
NAMES=['H_hi','H_ho','H_ci','H_co','mh','mc']
PEV=15.0; PB=5.0
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
H_ci=h3+(P('H','P',PEV*1e5,'S',s3,'R245fa')-h3)/0.75
H_hi=P('H','T',403.15,'P',PB*1e5,'Water')
def truth(three):
    mc=4.1321 if three else 4.7347
    H_co=P('H','T',395.15,'P',PEV*1e5,'R245fa') if three else P('H','P',PEV*1e5,'Q',1,'R245fa')
    H_ho=H_hi-mc*(H_co-H_ci)/10.0
    return dict(H_hi=H_hi,H_ho=H_ho,H_ci=H_ci,H_co=H_co,mh=10.0,mc=mc)
def run(unk,three):
    TR=truth(three)
    nw=Network(); nw.units.set_defaults(temperature='degC',pressure='bar',enthalpy='kJ/kg',mass_flow='kg/s',heat='kW')
    sb,kb,sw,kw=Source('sb'),Sink('kb'),Source('sw'),Sink('kw'); hx=MovingBoundaryHeatExchanger('evap')
    a=Connection(sb,'out1',hx,'in1'); b=Connection(hx,'out1',kb,'in1'); c=Connection(sw,'out1',hx,'in2'); d=Connection(hx,'out2',kw,'in1')
    nw.add_conns(a,b,c,d)
    attrs=dict(pr1=1,pr2=1,td_pinch=8.0)
    if three: attrs['ttd_u']=8.0
    hx.set_attr(**attrs)
    ca,cb,cc,cd={},{},{},{}
    if 'H_hi' not in unk: ca['T']=P('T','H',TR['H_hi'],'P',PB*1e5,'Water')-273.15
    if 'H_ho' not in unk: cb['T']=P('T','H',TR['H_ho'],'P',PB*1e5,'Water')-273.15
    if 'H_ci' not in unk: cc['h']=TR['H_ci']/1e3
    if 'H_co' not in unk:
        if three: cd['T']=P('T','H',TR['H_co'],'P',PEV*1e5,'R245fa')-273.15
        else: cd['x']=1
    if 'mh' not in unk: ca['m']=10.0
    if 'mc' not in unk: cc['m']=TR['mc']
    a.set_attr(fluid={'Water':1},p=PB,**ca); b.set_attr(**cb); c.set_attr(fluid={'R245fa':1},p=PEV,**cc); d.set_attr(**cd)
    t=time.perf_counter()
    try:
        nw.solve('design',init_only=False) if False else nw.solve('design')
        st=nw.status
        got=dict(H_hi=a.h.val*1e3,H_ho=b.h.val*1e3,H_ci=c.h.val*1e3,H_co=d.h.val*1e3,mh=a.m.val,mc=c.m.val)
        vals=[got[n] for n in unk]
        if st!=0 or any(v is None or not np.isfinite(v) for v in vals): res=('not solved',f'status {st}')
        else:
            e=max(abs(got[n]-TR[n])/abs(TR[n]) for n in unk)
            res=('solved' if e<2e-3 else 'solved, different answer',f'maxrel={e:.1e}')
    except Exception as ex: res=('raised',type(ex).__name__+': '+str(ex)[:80].replace('\n',' '))
    return res+(time.perf_counter()-t,)
out=[]
for k,three in ((2,False),(3,True)):
    for i,u in enumerate(itertools.combinations(NAMES,k),1):
        r=run(set(u),three); out.append(dict(k=k,i=i,u=list(u),status=r[0],detail=r[1],t=r[2])); print(k,i,','.join(u),r[0],r[1],f'{r[2]:.3f}s',flush=True)
json.dump(out,open('tespy_sweep.json','w'),indent=1)
