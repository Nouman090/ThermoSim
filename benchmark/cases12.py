import warnings, io, contextlib
warnings.filterwarnings('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
PEV=15e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
h4=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75
def solve(build, label):
    M=ThermodynamicModel(); hx=build(M)
    try:
        with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('always'); hx.Cal()
        pts=lambda k: M.Point[k]
        print(f'{label}\n  hot : in {pts("b_in").T-273.15:.3f} C  out {pts("b_out").T-273.15:.3f} C  m={pts("b_in").Mass_flowrate:.4f}'
              f'\n  cold: in {pts("w_in").T-273.15:.3f} C  out {pts("w_out").T-273.15:.3f} C (Q={pts("w_out").Q if hasattr(pts("w_out"),"Q") else "-"})  m={pts("w_in").Mass_flowrate:.4f}'
              f'\n  min_dT={hx.min_dT:.4f}  warnings={[str(x.message)[:90] for x in w]}')
    except Exception as e:
        print(f'{label}\n  EXCEPTION {type(e).__name__}: {str(e)[:260]}')
def case1(M):   # {H_hi,H_ho}: both flows known, cold side fully known
    M.add_point('Water','b_in',Mass_flowrate=10.0,P=5e5); M.add_point('Water','b_out',P=5e5)
    M.add_point('R245fa','w_in',Mass_flowrate=4.7347,P=PEV,H=h4); M.add_point('R245fa','w_out',P=PEV,Q=1)
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx; return hx
def case2(M):   # {H_ci,H_co}: hot side fully known, both flows known
    M.add_point('Water','b_in',Mass_flowrate=10.0,T=403.15,P=5e5); M.add_point('Water','b_out',T=104.78+273.15,P=5e5)
    M.add_point('R245fa','w_in',Mass_flowrate=4.7347,P=PEV); M.add_point('R245fa','w_out',P=PEV)
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx; return hx
solve(case1,'CASE 1  {H_hi,H_ho} unknown  (truth: hot 130.000 -> 104.78 C)')
solve(case2,'CASE 2  {H_ci,H_co} unknown  (truth: cold in = h4 (~48 C subcooled-side state), out = sat. vapour)')
print('truth h4 -> T =', P('T','H',h4,'P',PEV,'R245fa')-273.15, 'C ; T_sat(15bar)=',P('T','P',PEV,'Q',1,'R245fa')-273.15)
