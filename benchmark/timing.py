import time, statistics, warnings, io, contextlib, logging
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
from CoolProp.CoolProp import PropsSI as P
import ThermoSim as TS
from ThermoSim import ThermodynamicModel, HeatExchanger
from tespy.networks import Network
from tespy.components import Source, Sink, MovingBoundaryHeatExchanger
from tespy.connections import Connection
PEV=15e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
h4=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75

def ts(three):
    M=ThermodynamicModel()
    M.add_point('Water','b_in',Mass_flowrate=10.0,T=403.15,P=5e5); M.add_point('Water','b_out',P=5e5)
    M.add_point('R245fa','w_in',P=PEV,H=h4)
    M.add_point('R245fa','w_out',**(dict(P=PEV) if three else dict(P=PEV,Q=1)))
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx
    with contextlib.redirect_stdout(io.StringIO()): hx.Cal()
    return M.Point['w_in'].Mass_flowrate

def te(three):
    nw=Network(); nw.units.set_defaults(temperature='degC',pressure='bar',enthalpy='kJ/kg',mass_flow='kg/s',heat='kW'); nw.iterinfo=False
    sb,kb,sw,kw=Source('sb'),Sink('kb'),Source('sw'),Sink('kw'); hx=MovingBoundaryHeatExchanger('e')
    a=Connection(sb,'out1',hx,'in1'); b=Connection(hx,'out1',kb,'in1'); c=Connection(sw,'out1',hx,'in2'); d=Connection(hx,'out2',kw,'in1')
    nw.add_conns(a,b,c,d); hx.set_attr(pr1=1,pr2=1,td_pinch=8.0,**(dict(ttd_u=8.0) if three else {}))
    a.set_attr(fluid={'Water':1},T=130,p=5,m=10); c.set_attr(fluid={'R245fa':1},p=15,h=h4/1000)
    if not three: d.set_attr(x=1)
    with contextlib.redirect_stdout(io.StringIO()): nw.solve('design')
    return c.m.val

for name,f in (('ThermoSim',ts),('TESPy 0.11.3',te)):
    for three,lab in ((False,'P1 two-unknown'),(True,'P3 three-unknown')):
        f(three)  # warm-up
        t=[];
        for _ in range(7):
            s=time.perf_counter(); m=f(three); t.append(time.perf_counter()-s)
        print(f'{name:16} {lab:17} median {statistics.median(t):7.3f} s  (min {min(t):.3f}, max {max(t):.3f}, n=7)  m_wf={m:.4f}')
