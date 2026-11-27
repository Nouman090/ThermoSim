import warnings, io, contextlib
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
pev=15e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa'); h4=h3+(P('H','P',pev,'S',s3,'R245fa')-h3)/0.75
def ht(M): return {k:(p.H,p.T,p.P) for k,p in M.Point.items()}
def mf(M): return {k:p.Mass_flowrate for k,p in M.Point.items()}
def case(label,setup):
    M=ThermodynamicModel(); hx=setup(M); b=ht(M); bm=mf(M)
    try:
        with contextlib.redirect_stdout(io.StringIO()): hx.Cal(); r='SOLVED'; msg=''
    except Exception as e: r='refused'; msg=str(e)[:230].replace('\n',' ')
    a=ht(M); am=mf(M)
    ch=[k for k in b if b[k]!=a[k]]; chm=[k for k in bm if bm[k]!=am[k]]
    print(f'{label}: {r}; points whose H/T/P changed: {ch or "none"}; mass-flow entries changed: {chm or "none"}\n   {msg}')
def add(M,bo,wo,bf=10.0,wf=None):
    M.add_point('Water','b_in',Mass_flowrate=bf,T=403.15,P=5e5); M.add_point('Water','b_out',**bo)
    M.add_point('R245fa','w_in',Mass_flowrate=wf,P=pev,H=h4); M.add_point('R245fa','w_out',**wo)
def mk(M,ha=True):
    hx=HeatExchanger(M,'E',8.0,'Evaporator',ha,'b_in','b_out','w_in','w_out'); M.Component['E']=hx; return hx
case('A  specification crossing at the terminal (outlet 135 C)', lambda M:(add(M,dict(P=5e5),dict(P=pev,T=408.15)),mk(M))[1])
case('B  specification crossing at the terminal (outlet 128 C)', lambda M:(add(M,dict(P=5e5),dict(P=pev,T=401.15)),mk(M))[1])
case('C  both mass flows unknown (case 15)', lambda M:(add(M,dict(P=5e5,T=363.15),dict(P=pev,Q=1),bf=None),mk(M))[1])
case('D  knowns at hot end only (case 9)', lambda M:(M.add_point('Water','b_in',Mass_flowrate=10.0,T=403.15,P=5e5),M.add_point('Water','b_out',T=363.15,P=5e5),M.add_point('R245fa','w_in',P=pev),M.add_point('R245fa','w_out',P=pev,Q=1),mk(M))[-1])
case('E  feasible 2-unknown (control)', lambda M:(add(M,dict(P=5e5),dict(P=pev,Q=1)),mk(M))[1])
