import warnings, json
warnings.simplefilter('ignore')
from hrsg import plant
out={}
rows=[]
for T in (120,125,130,135,140,145,150,155,160):
    try:
        r=plant(60,T); rows.append((T,r['m'],r['mind'],r['eta'],r['Wnet']))
    except Exception as e: rows.append((T,'ERR',str(e)[:60]))
for x in rows: print(60,x)
# find stack T where min approach = 15 and = 0 (60 bar)
from scipy.optimize import brentq
f15=lambda T: plant(60,T)['mind']-15.0
f0=lambda T: plant(60,T)['mind']
T15=brentq(f15,120,200,xtol=1e-3); T0_=brentq(f0,120,200,xtol=1e-3)
r15=plant(60,T15); r0=plant(60,T0_)
print('60 bar: stack for 15K',round(T15,2),round(r15['m'],3),round(r15['eta'],3),round(r15['Wnet'],3))
print('60 bar: stack for 0K',round(T0_,2),round(r0['m'],3),round(r0['eta'],3))
# 40 bar B and best feasible across pressures with 15K design (stack chosen at min approach 15 K)
best=[]
for p in (20,30,40,50,60,70,80):
    try:
        Tp=brentq(lambda T: plant(p,T)['mind']-15.0,100,260,xtol=1e-3)
        r=plant(p,Tp); best.append((p,round(Tp,2),round(r['m'],3),round(r['eta'],3)))
    except Exception as e: best.append((p,'ERR',str(e)[:50]))
for b in best: print(b)
json.dump(dict(rows=rows,T15=T15,r15=[r15['m'],r15['eta'],r15['Wnet']],T0=T0_,best=best),open('feas.json','w'))
