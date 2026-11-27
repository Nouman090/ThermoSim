import warnings, logging, io, contextlib
warnings.filterwarnings('ignore')
exec(open('p234_tespy.py').read().split("# P2: terminals fixed")[0])
from ThermoSim import ThermodynamicModel, HeatExchanger
class H(logging.Handler):
    def __init__(s): super().__init__(logging.WARNING); s.rec=[]
    def emit(s,r): s.rec.append(r.getMessage())
for Tout in (100.0, 98.0):
    # ThermoSim
    M=ThermodynamicModel()
    M.add_point('Water','b_in',Mass_flowrate=10.0,T=403.15,P=5e5); M.add_point('Water','b_out',P=5e5,T=Tout+273.15)
    M.add_point('R245fa','w_in',P=15e5,H=h4*1000); M.add_point('R245fa','w_out',P=15e5,Q=1)
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx
    with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always'); hx.Cal()
    print(f'b_out={Tout} C | ThermoSim: min_dT={hx.min_dT:.2f} K, warnings={len(w)}', ('-> '+str(w[0].message)[:110]) if w else '')
    # TESPy
    h=H(); r=logging.getLogger(); r.addHandler(h)
    nw,hxt,a,b,c,d=build({}, dict(T=130,m=10), dict(T=Tout), dict(h=h4), dict(x=1)); nw.solve('design')
    r.removeHandler(h)
    print(f'            | TESPy:     td_pinch={hxt.td_pinch.val:.2f} K, status={nw.status}, warnings={len(h.rec)}')
