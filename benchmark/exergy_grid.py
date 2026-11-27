import warnings, json, sys
warnings.simplefilter('ignore')
from orc import run
def metrics(pev, ppt):
    M=run(pev=pev, ppt=ppt); ev=M.Component['Evaporator']
    b1,b2=M.Point['b1'],M.Point['b2']
    dEx=(b1.Ex-b2.Ex)
    return dict(p=pev,ppt=ppt,m=M.Point['1'].Mass_flowrate,Tb2=b2.T-273.15,W=M.Net_power/1e3,eta=M.Efficiency,
        ExD={k:c.Ex_D/1e3 for k,c in M.Component.items()},dEx=dEx/1e3,psi=100*M.Net_power/dEx,mindT=ev.min_dT,Wt=M.Power_Out/1e3,Wp=M.Power_In/1e3)
res={}
res['base']=metrics(15,8)
res['ppt_sweep']=[metrics(15,a) for a in (4,6,8,10,12)]
grid=[]
for p in (8,10,12,15,18):
    for a in (4,6,8,10,12):
        try: grid.append(metrics(p,a))
        except Exception as e: grid.append(dict(p=p,ppt=a,err=str(e)[:80]))
res['grid']=grid
# feasibility limit
lim={}
for a in (4,8,12):
    lo,hi=8.0,30.0
    ok=lambda p:(lambda r:r)(True)
    def feas(p):
        try:
            m=metrics(p,a); return m['mindT']>=a-0.05 and m['m']>0
        except Exception: return False
    if not feas(lo): lim[a]=None; continue
    if feas(hi): lim[a]=hi; continue
    for _ in range(12):
        mid=(lo+hi)/2
        if feas(mid): lo=mid
        else: hi=mid
    lim[a]=round(lo,2)
res['limit']=lim
json.dump(res,open('exergy_grid.json','w'),indent=1)
print(json.dumps(res['base'],indent=1)); 
for r in res['ppt_sweep']: print(r['ppt'],round(r['m'],3),round(r['Tb2'],2),round(r['W'],1),round(r['eta'],2),{k:round(v,1) for k,v in r['ExD'].items()},round(r['psi'],1))
print(lim)
for r in grid: print(r.get('p'),r.get('ppt'),r.get('err') or (round(r['m'],3),round(r['W'],1),round(r['eta'],2),round(r['psi'],1),round(r['mindT'],2)))
