import warnings, json
warnings.simplefilter('ignore')
import numpy as np
from CoolProp.CoolProp import PropsSI as P
from hrsg import plant, G
def naive(r, p_bar, N, hot='Air', hotP=1e5):
    M=r['M']; hx=r['hx']
    gi,go,wi,wo=M.Point['g_in'],M.Point['g_out'],M.Point['w_in'],M.Point['w_out']
    m_g=100.0; m_w=M.Point['w_in'].Mass_flowrate
    Q=m_g*(gi.H-go.H)
    mins=1e9
    for k in range(N+1):
        q=Q*k/N
        Hc=wi.H+q/m_w; Hh=go.H+q/m_g
        Tc=P('T','H',Hc,'P',p_bar*1e5,'Water'); Th=P('T','H',Hh,'P',hotP,'Air')
        mins=min(mins,Th-Tc)
    return mins
res={'A':[], 'B':[]}
for name,(p,T) in dict(A=(60,120),B=(40,150)).items():
    for N in (5,10,20,50,100,200,500,1000):
        r=plant(p,T,N=N)
        res[name].append((N,r['mind'],naive(r,p,N),r['m']))
        print(name,N,round(r['mind'],4),round(res[name][-1][2],4),round(r['m'],5),flush=True)
json.dump(res,open('disc.json','w'))
