import warnings, io, contextlib, json, time, numpy as np
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
# CO2 gas cooler: hot = CO2 at 100 bar, 120 C -> 32 C, 1 kg/s; cold = water at 3 bar entering 15 C; cold outlet and cold flow unknown; DT_pp = 5 K
PH=100e5; PC=3e5; mh=1.0; PPT=5.0
def solve(N):
    M=ThermodynamicModel()
    M.add_point('CO2','hi',Mass_flowrate=mh,T=393.15,P=PH); M.add_point('CO2','ho',T=305.15,P=PH)
    M.add_point('Water','ci',T=288.15,P=PC); M.add_point('Water','co',P=PC)
    hx=HeatExchanger(M,'G',PPT,'double_pipe',False,'hi','ho','ci','co',div_N=N); M.Component['G']=hx
    t=time.perf_counter()
    with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always'); hx.Cal()
    return M,hx,time.perf_counter()-t,len(w)
def fine(M,K=40001):
    hi,ho,ci,co=[M.Point[k] for k in ('hi','ho','ci','co')]
    mc=M.Point['ci'].Mass_flowrate; Q=mh*(hi.H-ho.H)
    a=np.linspace(0,1,K); Th=np.array([P('T','H',ho.H+x*Q/mh,'P',PH,'CO2') for x in a]); Tc=np.array([P('T','H',ci.H+x*Q/mc,'P',PC,'Water') for x in a])
    d=Th-Tc; i=int(d.argmin()); return d[i],a[i],Th[i]-273.15
rows=[]
for N in (5,10,20,50,100,200,500):
    try:
        M,hx,dt,nw=solve(N); mc=M.Point['ci'].Mass_flowrate; Tco=M.Point['co'].T-273.15
        dmin,alpha,Tpin=fine(M)
        rows.append(dict(N=N,mc=mc,Tco=Tco,fine_min=dmin,alpha=alpha,Tpin=Tpin,t=dt,warns=nw)); print(N,round(mc,5),round(Tco,3),round(dmin,4),round(alpha,3),round(Tpin,1),round(dt,1),nw,flush=True)
    except Exception as e: rows.append(dict(N=N,err=str(e)[:100])); print(N,'ERR',str(e)[:160],flush=True)
json.dump(rows,open('conv.json','w'),indent=1)
