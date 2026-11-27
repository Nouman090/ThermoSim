"""Sweep every 2- and 3-unknown combination of {H_hi,H_ho,H_ci,H_co,mh,mc} on the ORC evaporator
(R245fa, brine 130 C / 10 kg/s, 8 K pinch, 15 bar) and record solved / refused / wrong."""
import warnings, io, contextlib, itertools, signal, time, sys, json
warnings.filterwarnings('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
PEV=15e5; PB=5e5; PPT=8.0
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
H_ci=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75
H_co=P('H','P',PEV,'Q',1,'R245fa'); mc=4.7347
H_hi=P('H','T',403.15,'P',PB,'Water'); mh=10.0
Q=mc*(H_co-H_ci); H_ho=H_hi-Q/mh          # exactly energy-consistent truth
TRUTH=dict(H_hi=H_hi,H_ho=H_ho,H_ci=H_ci,H_co=H_co,mh=mh,mc=mc)
NAMES=['H_hi','H_ho','H_ci','H_co','mh','mc']
class TO(Exception): pass
def _alarm(*a): raise TO()
signal.signal(signal.SIGALRM,_alarm)
def run(unk):
    M=ThermodynamicModel()
    kw=lambda n,fl,P_: ({} if n not in TRUTH or n in unk else {})
    def pt(name,fluid,Pr,hname,mname):
        d=dict(P=Pr)
        if hname not in unk: d['H']=TRUTH[hname]
        if mname not in unk: d['Mass_flowrate']=TRUTH[mname]
        M.add_point(fluid,name,**d)
    pt('hi','Water',PB,'H_hi','mh'); pt('ho','Water',PB,'H_ho','mh')
    pt('ci','R245fa',PEV,'H_ci','mc'); pt('co','R245fa',PEV,'H_co','mc')
    hx=HeatExchanger(M,'E',PPT,'Evaporator',False,'hi','ho','ci','co'); M.Component['E']=hx
    t=time.perf_counter(); signal.alarm(150)
    try:
        with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('always'); hx.Cal()
        signal.alarm(0)
        got=dict(H_hi=M.Point['hi'].H,H_ho=M.Point['ho'].H,H_ci=M.Point['ci'].H,H_co=M.Point['co'].H,
                 mh=M.Point['hi'].Mass_flowrate,mc=M.Point['ci'].Mass_flowrate)
        err={n:abs(got[n]-TRUTH[n])/abs(TRUTH[n]) for n in unk}
        return dict(status='solved',time=time.perf_counter()-t,min_dT=hx.min_dT,maxrel=max(err.values()),
                    warn=len(w),got={n:got[n] for n in unk})
    except TO: return dict(status='timeout',time=150)
    except Exception as e:
        signal.alarm(0); return dict(status='refused',time=time.perf_counter()-t,msg=str(e)[:170].replace('\n',' '))
k=int(sys.argv[1]); out=[]
for i,unk in enumerate(itertools.combinations(NAMES,k),1):
    r=run(set(unk)); r['case']=i; r['unknowns']=list(unk); out.append(r)
    tail=(f"maxrel={r['maxrel']:.2e} minDT={r['min_dT']:.3f} warn={r['warn']}" if r['status']=='solved' else r.get('msg',''))
    print(f"{k}-unk #{i:2d} {{{','.join(unk)}}}: {r['status']:8s} {r['time']:5.1f}s {tail}",flush=True)
json.dump(out,open(f'sweep_{k}.json','w'),indent=1,default=float)
