import json, warnings; warnings.simplefilter('ignore')
from chena_s import run
cases=[('base',{}),
('P_evap +0.1 bar',dict(dP=0.1)),('P_evap -0.1 bar',dict(dP=-0.1)),
('T_geo,in +0.5 K',dict(Tgi=73.83)),('T_geo,in -0.5 K',dict(Tgi=72.83)),
('T_geo,out +0.5 K',dict(Tgo=54.94)),('T_geo,out -0.5 K',dict(Tgo=53.94)),
('kappa 1.03',dict(kap=1.03)),('kappa 1.05',dict(kap=1.05)),
('T_4 +0.5 K',dict(Te_in=13.35)),('T_4 -0.5 K',dict(Te_in=12.35)),
('P_water 1 bar',dict(Pw=1e5)),('P_water 5 bar',dict(Pw=5e5)),
('P_cond 4.5 bar',dict(Pc=4.5)),('eta_t 0.85',dict(eta=0.85))]
res={}
for Pev in (16.0,16.95):
    for nm,kw in cases:
        kw=dict(kw); dP=kw.pop('dP',0.0)
        try: r=run(Pev=Pev+dP,N=100,**kw)
        except Exception as e: r={'err':str(e)[:80]}
        res[f'{Pev}|{nm}']=r
        print(Pev,nm,{k:(round(float(v),4) if not isinstance(v,str) else v) for k,v in r.items()},flush=True)
json.dump(res,open('chena_sens.json','w'),indent=1,default=float)
