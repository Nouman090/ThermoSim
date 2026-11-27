import warnings, io, contextlib, json, numpy as np
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
out={}
# --- C3: same-side example (water as near-constant-cp fluid, 1 kg/s each, cold 20->60 C, PPT 10 K), unknown H_hi,H_ho
M=ThermodynamicModel()
M.add_point('Water','hi',P=5e5,Mass_flowrate=1.0); M.add_point('Water','ho',P=5e5)
M.add_point('Water','ci',P=5e5,T=293.15,Mass_flowrate=1.0); M.add_point('Water','co',P=5e5,T=333.15)
hx=HeatExchanger(M,'X',10.0,'double_pipe',None,'hi','ho','ci','co'); M.Component['X']=hx
with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
    warnings.simplefilter('always'); hx.Cal()
out['c3']=dict(Thi=M.Point['hi'].T-273.15,Tho=M.Point['ho'].T-273.15,mind=hx.min_dT,warn=len(w),cp20=P('C','T',293.15,'P',5e5,'Water'),cp60=P('C','T',333.15,'P',5e5,'Water'))
print(out['c3'])
# --- C4: family of feasible designs, ORC evaporator, hot-end outlet temperature T_co swept (2-unknown solves), plus anchor (3-unknown)
PEV=15e5;PB=5e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa'); H_ci=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75
H_hi=P('H','T',403.15,'P',PB,'Water')
def solve(Tco=None):
    M=ThermodynamicModel()
    M.add_point('Water','hi',P=PB,T=403.15,Mass_flowrate=10.0); M.add_point('Water','ho',P=PB)
    M.add_point('R245fa','ci',P=PEV,H=H_ci)
    M.add_point('R245fa','co',P=PEV,**({'Q':1} if Tco is None else {'T':Tco+273.15}))
    hx=HeatExchanger(M,'X',8.0,'Evaporator',False,'hi','ho','ci','co',div_N=40); M.Component['X']=hx
    with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always'); hx.Cal()
    return M.Point['ci'].Mass_flowrate, M.Point['ho'].T-273.15, hx.min_dT, M.Point['hi'].T-M.Point['co'].T
Tsat=P('T','P',PEV,'Q',1,'R245fa')-273.15
fam=[]
for T in (None,110,114,118,121):
    m,Tho,md,dhot=solve(T); fam.append(dict(Tco=(Tsat if T is None else T),m=m,Tho=Tho,mind=md,dhot=dhot)); print(fam[-1])
out['family']=fam
# --- C4: Jacobian of the three closure equations at the 3-unknown solution, unknowns x=(mc,H_co,H_ho)
mh=10.0; PPT=8.0
Hc_s=P('H','P',PEV,'Q',0,'R245fa')
def resid(x,N=400):
    mc,Hco,Hho=x
    Th=lambda h:P('T','H',h,'P',PB,'Water'); Tc=lambda h:P('T','H',h,'P',PEV,'R245fa')
    g1=mh*(H_hi-Hho)-mc*(Hco-H_ci)
    a=np.linspace(0,1,N+1); hs=Hho+a*(H_hi-Hho); cs=H_ci+a*(Hco-H_ci)
    # insert exact bubble point
    ab=(Hc_s-H_ci)/(Hco-H_ci)
    a=np.sort(np.append(a,ab)); hs=Hho+a*(H_hi-Hho); cs=H_ci+a*(Hco-H_ci)
    g2=min(Th(h)-Tc(c) for h,c in zip(hs,cs))-PPT
    g3=Th(H_hi)-Tc(Hco)-PPT
    return np.array([g1/1e5,g2,g3])
x0=np.array([4.1321,P('H','T',395.15,'P',PEV,'R245fa'),H_hi-4.1321*(P('H','T',395.15,'P',PEV,'R245fa')-H_ci)/mh])
scale=np.array([1e-3,1e2,1e2]); J=np.zeros((3,3))
for j in range(3):
    d=np.zeros(3); d[j]=scale[j]*1.0
    J[:,j]=(resid(x0+d)-resid(x0-d))/(2*d[j])
sv=np.linalg.svd(J*np.array([1.0,1e5,1e5]),compute_uv=False)  # unit-free-ish scaling by typical magnitudes
Jn=J*x0  # column scaling by variable magnitude
sv2=np.linalg.svd(Jn,compute_uv=False)
print('residual at x0',resid(x0)); print('J\n',J); print('sv (scaled by x)',sv2,'cond',sv2[0]/sv2[-1],'det',np.linalg.det(J))
out['jac']=dict(res=resid(x0).tolist(),J=J.tolist(),sv=sv2.tolist(),cond=float(sv2[0]/sv2[-1]))
json.dump(out,open('c345.json','w'),indent=1,default=float)
