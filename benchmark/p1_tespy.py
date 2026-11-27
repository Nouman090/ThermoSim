import time, warnings
from CoolProp.CoolProp import PropsSI as P
from tespy.networks import Network
from tespy.components import Source, Sink, MovingBoundaryHeatExchanger
from tespy.connections import Connection

pev=15.0
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
h4=h3+(P('H','P',pev*1e5,'S',s3,'R245fa')-h3)/0.75

nw=Network()
nw.units.set_defaults(temperature='degC', pressure='bar', enthalpy='kJ/kg', mass_flow='kg/s', heat='kW')
nw.set_attr(iterinfo=False)
src_b=Source('b_src'); snk_b=Sink('b_snk'); src_w=Source('w_src'); snk_w=Sink('w_snk')
hx=MovingBoundaryHeatExchanger('evap')
c1=Connection(src_b,'out1',hx,'in1',label='b_in')
c2=Connection(hx,'out1',snk_b,'in1',label='b_out')
c3=Connection(src_w,'out1',hx,'in2',label='w_in')
c4=Connection(hx,'out2',snk_w,'in1',label='w_out')
nw.add_conns(c1,c2,c3,c4)
hx.set_attr(pr1=1,pr2=1,td_pinch=8.0)
c1.set_attr(fluid={'Water':1},T=130,p=5,m=10)
c3.set_attr(fluid={'R245fa':1},p=pev,h=h4/1000)
c4.set_attr(x=1)
t=time.perf_counter()
try:
    nw.solve('design')
    dt=time.perf_counter()-t
    print('status',nw.status, 'time',round(dt,3),'s')
    print('m_wf =',c3.m.val,'kg/s | T_b_out =',c2.T.val,'C | Q =',hx.Q.val,'kW | td_pinch =',hx.td_pinch.val)
except Exception as e:
    print('FAIL',type(e).__name__,str(e)[:300])
