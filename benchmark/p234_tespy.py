import time, logging, io, warnings
warnings.filterwarnings('ignore')
from CoolProp.CoolProp import PropsSI as P
from tespy.networks import Network
from tespy.components import Source, Sink, MovingBoundaryHeatExchanger
from tespy.connections import Connection

pev=15.0
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
h4=(h3+(P('H','P',pev*1e5,'S',s3,'R245fa')-h3)/0.75)/1000

def build(hx_attrs, c1=None, c2=None, c3=None, c4=None):
    nw=Network(); nw.units.set_defaults(temperature='degC',pressure='bar',enthalpy='kJ/kg',mass_flow='kg/s',heat='kW')
    sb,kb,sw,kw=Source('sb'),Sink('kb'),Source('sw'),Sink('kw')
    hx=MovingBoundaryHeatExchanger('evap')
    a=Connection(sb,'out1',hx,'in1',label='b_in'); b=Connection(hx,'out1',kb,'in1',label='b_out')
    c=Connection(sw,'out1',hx,'in2',label='w_in'); d=Connection(hx,'out2',kw,'in1',label='w_out')
    nw.add_conns(a,b,c,d)
    hx.set_attr(pr1=1,pr2=1,**hx_attrs)
    a.set_attr(fluid={'Water':1},p=5,**(c1 or {}))
    b.set_attr(**(c2 or {}))
    c.set_attr(fluid={'R245fa':1},p=pev,**(c3 or {}))
    d.set_attr(**(c4 or {}))
    return nw,hx,a,b,c,d

def run(label, *args, **kw):
    log=io.StringIO(); h=logging.StreamHandler(log); h.setLevel(logging.DEBUG)
    lg=logging.getLogger('root'); lg.addHandler(h)
    for n in ('tespy',): logging.getLogger(n).addHandler(h)
    t=time.perf_counter(); status=''
    try:
        nw,hx,a,b,c,d=build(*args,**kw)
        nw.solve('design'); status=f'solver status={nw.status}'
        res=(f'm_wf={c.m.val:.4f}  T_b_in={a.T.val:.2f}  T_b_out={b.T.val:.2f}  T_w_out={d.T.val:.2f}  '
             f'm_b={a.m.val:.3f}  td_pinch={getattr(hx.td_pinch,"val",None)}')
    except Exception as e:
        res=f'EXCEPTION {type(e).__name__}: {str(e)[:350]}'; status='raised'
    dt=time.perf_counter()-t
    print(f'--- {label} [{status}] {dt:.2f}s\n   {res}')
    msgs=[l for l in log.getvalue().splitlines() if l.strip()]
    for l in msgs[:6]: print('   log:',l[:300])

# P2: terminals fixed, no pinch spec (brine out 90 C) -> does TESPy flag the crossing?
run('P2 crossing: b_out=90C, no td_pinch spec', {}, dict(T=130,m=10), dict(T=90), dict(h=h4), dict(x=1))
# P2b: same but with td_pinch=8 specified as a constraint (over-specified: m_wf now fixed by EB AND pinch)
run('P2b crossing: b_out=90C AND td_pinch=8', dict(td_pinch=8.0), dict(T=130,m=10), dict(T=90), dict(h=h4), dict(x=1))
# P3: three unknowns {b_out, w_out T, m_wf}: td_pinch=8 + upper terminal difference 8
run('P3 three unknowns: td_pinch=8 + ttd_u=8', dict(td_pinch=8.0, ttd_u=8.0), dict(T=130,m=10), {}, dict(h=h4), {})
# P4a: both flows unknown
run('P4a unsupported: both flows unknown', dict(td_pinch=8.0), dict(T=130), dict(T=90), dict(h=h4), dict(x=1))
# P4b: hot side T unknown at both ends, flows known, cold side known
run('P4b {H_hi,H_ho} unknown: td_pinch=8, flows known', dict(td_pinch=8.0), dict(m=10), {}, dict(h=h4,m=4.7347), dict(x=1))
