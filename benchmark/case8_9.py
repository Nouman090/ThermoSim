import warnings, io, contextlib, logging
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
exec(open('p234_tespy.py').read().split("# P2: terminals fixed")[0])   # build(), h4 [kJ/kg], pev
PEV=15e5; h4J=h4*1000
# ---- ThermoSim: case 8 = {H_hi,H_co} unknown; knowns H_ho (cold end), H_ci (cold end), mh, mc
M=ThermodynamicModel()
M.add_point('Water','b_in',Mass_flowrate=10.0,P=5e5); M.add_point('Water','b_out',T=104.78+273.15,P=5e5)
M.add_point('R245fa','w_in',Mass_flowrate=4.7347,P=PEV,H=h4J); M.add_point('R245fa','w_out',P=PEV)
hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx
try:
    with contextlib.redirect_stdout(io.StringIO()): hx.Cal()
    print('ThermoSim case 8: SOLVED', M.Point['b_in'].T-273.15, M.Point['w_out'].T-273.15, hx.min_dT)
except Exception as e: print('ThermoSim case 8: EXCEPTION', str(e)[:300])
# ---- ThermoSim: case 9 = {H_ho,H_ci} unknown; knowns H_hi (hot end), H_co (hot end), mh, mc
M=ThermodynamicModel()
M.add_point('Water','b_in',Mass_flowrate=10.0,T=403.15,P=5e5); M.add_point('Water','b_out',P=5e5)
M.add_point('R245fa','w_in',Mass_flowrate=4.7347,P=PEV); M.add_point('R245fa','w_out',P=PEV,Q=1)
hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx
try:
    with contextlib.redirect_stdout(io.StringIO()): hx.Cal()
    print('ThermoSim case 9: SOLVED', M.Point['b_out'].T-273.15, M.Point['w_in'].T-273.15, hx.min_dT)
except Exception as e: print('ThermoSim case 9: EXCEPTION', str(e)[:300])
# ---- TESPy case 8: hot out T, cold in h, both m known, td_pinch=8 ; hot-in & cold-out free
nw,hxt,a,b,c,d=build(dict(td_pinch=8.0), dict(m=10), dict(T=104.78), dict(h=h4,m=4.7347), {})
try: nw.solve('design'); print(f'TESPy case 8: status={nw.status}  T_hot_in={a.T.val:.3f}  T_cold_out={d.T.val:.3f}  pinch={hxt.td_pinch.val:.3f}')
except Exception as e: print('TESPy case 8: EXCEPTION', str(e)[:250])
# ---- TESPy case 9: hot in T, cold out x=1, both m known, td_pinch=8 ; hot-out & cold-in free
nw,hxt,a,b,c,d=build(dict(td_pinch=8.0), dict(T=130,m=10), {}, dict(m=4.7347), dict(x=1))
try: nw.solve('design'); print(f'TESPy case 9: status={nw.status}  T_hot_out={b.T.val:.3f}  T_cold_in={c.T.val:.3f}  pinch={hxt.td_pinch.val:.3f}')
except Exception as e: print('TESPy case 9: EXCEPTION', str(e)[:250])
