import warnings, io, contextlib, signal, time, json
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
def a(*x): raise TimeoutError
signal.signal(signal.SIGALRM,a)
def h4(pev):
    h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
    return h3+(P('H','P',pev,'S',s3,'R245fa')-h3)/0.75
def three(Tb, pev_bar):
    pev=pev_bar*1e5; M=ThermodynamicModel()
    M.add_point('Water','b_in',Mass_flowrate=10.0,T=Tb+273.15,P=5e5); M.add_point('Water','b_out',P=5e5)
    M.add_point('R245fa','w_in',P=pev,H=h4(pev)); M.add_point('R245fa','w_out',P=pev)
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx
    t=time.perf_counter(); signal.alarm(120)
    with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always'); hx.Cal()
    signal.alarm(0)
    m=M.Point['w_in'].Mass_flowrate; bi,bo,wi,wo=[M.Point[k] for k in ('b_in','b_out','w_in','w_out')]
    # independent evaluation of the two approaches from the returned states
    d_hot=bi.T-wo.T
    Tsat=P('T','P',pev,'Q',0,'R245fa'); hl=P('H','P',pev,'Q',0,'R245fa')
    Qb=m*(hl-wi.H); Tb_onset=P('T','H',bo.H+Qb/10.0,'P',5e5,'Water')
    d_int=Tb_onset-Tsat
    return dict(Tb=Tb,p=pev_bar,m=m,Tbo=bo.T-273.15,Two=wo.T-273.15,d_hot=d_hot,d_int=d_int,mind=hx.min_dT,t=time.perf_counter()-t,warns=len(w))
rows=[]
for Tb,p in ((130,15),(130,12),(130,18),(120,12),(120,15),(140,15),(140,20),(150,20)):
    try: r=three(Tb,p); rows.append(r); print({k:round(v,4) for k,v in r.items()},flush=True)
    except Exception as e: signal.alarm(0); rows.append(dict(Tb=Tb,p=p,err=str(e)[:80])); print(Tb,p,'ERR',str(e)[:100],flush=True)
json.dump(rows,open('dp.json','w'))
