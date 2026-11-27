import warnings, io, contextlib, json, numpy as np
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
from scipy.optimize import brentq
FR=(0.006,0.02,0.05,0.1,0.2,0.3)
# (a) ORC evaporator, three unknowns {H_ho,H_co,mc}
pev=15e5; h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa'); h4=h3+(P('H','P',pev,'S',s3,'R245fa')-h3)/0.75
resa=[]
for fr in FR:
    M=ThermodynamicModel()
    M.add_point('Water','b_in',Mass_flowrate=10.0,T=403.15,P=5e5); M.add_point('Water','b_out',P=5e5)
    M.add_point('R245fa','w_in',P=pev,H=h4); M.add_point('R245fa','w_out',P=pev)
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out',anchor_exclude_frac=fr); M.Component['E']=hx
    with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always'); hx.Cal()
    m=M.Point['w_in'].Mass_flowrate; bo=M.Point['b_out']; wo=M.Point['w_out']
    resa.append(dict(frac=fr,m=m,Tbo=bo.T-273.15,Two=wo.T-273.15,mind=hx.min_dT,warns=len(w)))
    print('ORC',fr,round(m,5),round(bo.T-273.15,3),round(wo.T-273.15,3),round(hx.min_dT,4),len(w),flush=True)
# alpha location of the boiling onset in the 0.1-band question
M0=None
# (b) condenser truth case
PH=3e5; PC=2e5; T_hi=343.15; T_ho=308.15; mh=2.0; T_ci=T_ho-8.0
H_hi=P('H','T',T_hi,'P',PH,'R245fa'); H_ho=P('H','T',T_ho,'P',PH,'R245fa'); H_ci=P('H','T',T_ci,'P',PC,'Water')
Q=mh*(H_hi-H_ho); H_dew=P('H','P',PH,'Q',1,'R245fa'); T_dew=P('T','P',PH,'Q',1,'R245fa'); Q_dew=mh*(H_dew-H_ho)
mc=brentq(lambda m: T_dew-P('T','H',H_ci+Q_dew/m,'P',PC,'Water')-8.0,5,200); H_co=H_ci+Q/mc
TR=dict(H_hi=H_hi,H_ho=H_ho,H_ci=H_ci,H_co=H_co,mh=mh,mc=mc); alpha_dew=Q_dew/Q
print('condenser dew-point alpha from cold end =',round(alpha_dew,4))
resb=[]
unk={'H_ci','H_hi','mh'}
for fr in FR:
    M=ThermodynamicModel()
    def pt(n,f_,Pr,h,m):
        d=dict(P=Pr)
        if h not in unk: d['H']=TR[h]
        if m not in unk: d['Mass_flowrate']=TR[m]
        M.add_point(f_,n,**d)
    pt('hi','R245fa',PH,'H_hi','mh'); pt('ho','R245fa',PH,'H_ho','mh'); pt('ci','Water',PC,'H_ci','mc'); pt('co','Water',PC,'H_co','mc')
    hx=HeatExchanger(M,'C',8.0,'Condenser',None,'hi','ho','ci','co',anchor_exclude_frac=fr); M.Component['C']=hx
    try:
        with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('always'); hx.Cal()
        got=dict(H_hi=M.Point['hi'].H,H_ci=M.Point['ci'].H,mh=M.Point['hi'].Mass_flowrate)
        e=max(abs(got[n]-TR[n])/abs(TR[n]) for n in unk)
        resb.append(dict(frac=fr,err=e,mind=hx.min_dT,warns=len(w))); print('COND',fr,f'{e:.1e}',round(hx.min_dT,4),len(w),flush=True)
    except Exception as ex:
        resb.append(dict(frac=fr,err=None,msg=str(ex)[:100])); print('COND',fr,'REFUSED',str(ex)[:100],flush=True)
json.dump(dict(orc=resa,cond=resb,alpha_dew=alpha_dew),open('thr.json','w'),indent=1)
