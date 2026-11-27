import warnings, io, contextlib, json
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from hrsg import G
Tin=G['T4']
def three(p_bar, ppt=15.0):
    M=ThermodynamicModel()
    hs=P('H','P',0.008e6,'Q',0,'Water'); ss=P('S','P',0.008e6,'Q',0,'Water'); hp=hs+(P('H','P',p_bar*1e5,'S',ss,'Water')-hs)/0.80
    M.add_point('Air','g_in',P=1e5,T=Tin,Mass_flowrate=100.0); M.add_point('Air','g_out',P=1e5)
    M.add_point('Water','w_in',P=p_bar*1e5,H=hp); M.add_point('Water','w_out',P=p_bar*1e5)
    hx=HeatExchanger(M,'HRSG',ppt,'Evaporator',True,'g_in','g_out','w_in','w_out'); M.Component['HRSG']=hx
    with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always'); hx.Cal()
    gi,go,wi,wo=[M.Point[k] for k in ('g_in','g_out','w_in','w_out')]
    m=wi.Mass_flowrate
    # independent two-pinch check: hot end, and bubble point
    d_hot=gi.T-wo.T
    hl=P('H','P',p_bar*1e5,'Q',0,'Water'); Tsat=P('T','P',p_bar*1e5,'Q',0,'Water')
    Qb=m*(hl-wi.H); Tg_b=P('T','H',go.H+Qb/100.0,'P',1e5,'Air'); d_bub=Tg_b-Tsat
    h5=wo.H; s5=P('S','H',h5,'P',p_bar*1e5,'Water'); h6=h5-0.88*(h5-P('H','P',0.008e6,'S',s5,'Water'))
    Wst=m*(h5-h6); Wp=m*(hp-hs); Wnet=G['Wt']-G['Wc']+Wst-Wp
    return dict(p=p_bar,Tco=wo.T-273.15,Tstack=go.T-273.15,m=m,d_hot=d_hot,d_bub=d_bub,mind=hx.min_dT,eta=100*Wnet/G['Q'],Wnet=Wnet/1e6,warns=len(w))
rows=[three(p) for p in (40,60)]
for r in rows: print({k:round(v,4) for k,v in r.items()})
json.dump(rows,open('hrsg3.json','w'),indent=1)
