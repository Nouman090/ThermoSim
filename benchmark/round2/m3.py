import warnings, io, contextlib, itertools, time, json, sys
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from scipy.optimize import brentq
from ThermoSim import ThermodynamicModel, HeatExchanger
NAMES=['H_hi','H_ho','H_ci','H_co','mh','mc']
# evaporator truth (hot-end anchor, 8 K)
PEV=15e5; PB=5e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
H_ci=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75
H_hi=P('H','T',403.15,'P',PB,'Water'); H_co=P('H','T',395.15,'P',PEV,'R245fa'); mh=10.0; mc=4.1321
EV=dict(truth=dict(H_hi=H_hi,H_ho=H_hi-mc*(H_co-H_ci)/mh,H_ci=H_ci,H_co=H_co,mh=mh,mc=mc),
        fl=('Water',PB,'R245fa',PEV),typ='Evaporator',HA=False)
# condenser truth (cold-end anchor + dew point pinch)
PH=3e5; PC=2e5; T_hi=343.15; T_ho=308.15; mh2=2.0; T_ci=T_ho-8.0
cH_hi=P('H','T',T_hi,'P',PH,'R245fa'); cH_ho=P('H','T',T_ho,'P',PH,'R245fa'); cH_ci=P('H','T',T_ci,'P',PC,'Water')
Q=mh2*(cH_hi-cH_ho); H_dew=P('H','P',PH,'Q',1,'R245fa'); T_dew=P('T','P',PH,'Q',1,'R245fa'); Qd=mh2*(H_dew-cH_ho)
mc2=brentq(lambda m: T_dew-P('T','H',cH_ci+Qd/m,'P',PC,'Water')-8.0,5,200)
CO=dict(truth=dict(H_hi=cH_hi,H_ho=cH_ho,H_ci=cH_ci,H_co=cH_ci+Q/mc2,mh=mh2,mc=mc2),fl=('R245fa',PH,'Water',PC),typ='Condenser',HA=None)
def run(spec,unk,N=20):
    TR=spec['truth']; fh,Ph,fc,Pc=spec['fl']; M=ThermodynamicModel()
    def pt(n,f,Pr,h,m):
        d=dict(P=Pr)
        if h not in unk: d['H']=TR[h]
        if m not in unk: d['Mass_flowrate']=TR[m]
        M.add_point(f,n,**d)
    pt('hi',fh,Ph,'H_hi','mh'); pt('ho',fh,Ph,'H_ho','mh'); pt('ci',fc,Pc,'H_ci','mc'); pt('co',fc,Pc,'H_co','mc')
    hx=HeatExchanger(M,'X',8.0,spec['typ'],spec['HA'],'hi','ho','ci','co',div_N=N); M.Component['X']=hx
    t=time.perf_counter()
    try:
        with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('always'); hx.Cal()
        got=dict(H_hi=M.Point['hi'].H,H_ho=M.Point['ho'].H,H_ci=M.Point['ci'].H,H_co=M.Point['co'].H,mh=M.Point['hi'].Mass_flowrate,mc=M.Point['ci'].Mass_flowrate)
        e=max(abs(got[n]-TR[n])/abs(TR[n]) for n in unk)
        return ('solved',time.perf_counter()-t,f'maxrel={e:.1e} minDT={hx.min_dT:.3f} warns={len(w)}')
    except Exception as ex: return ('refused',time.perf_counter()-t,type(ex).__name__+': '+str(ex)[:110].replace('\n',' '))
if __name__=='__main__':
    out=[]
    for spec in (EV,CO):
        for i,unk in enumerate(itertools.combinations(NAMES,3),1):
            r=run(spec,set(unk)); out.append((spec['typ'],i,list(unk),r)); print(spec['typ'][:4],i,','.join(unk),r[0],f'{r[1]:.1f}s',r[2],flush=True)
    json.dump(out,open('m3.json','w'),default=float)
