import warnings, io, contextlib, copy, time
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
PEV=15e5; PB=5e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
H_ci=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75
H_co=P('H','P',PEV,'Q',1,'R245fa'); mc=4.7347; H_hi=P('H','T',403.15,'P',PB,'Water'); mh=10.0
H_ho=H_hi-mc*(H_co-H_ci)/mh
TR=dict(H_hi=H_hi,H_ho=H_ho,H_ci=H_ci,H_co=H_co,mh=mh,mc=mc)
def snap(M): return {k:{a:v for a,v in vars(p).items() if isinstance(v,(int,float,str,type(None)))} for k,p in M.Point.items()}
def trial(unk,N=20,label=''):
    M=ThermodynamicModel()
    def pt(n,f,Pr,h,m):
        d=dict(P=Pr)
        if h not in unk: d['H']=TR[h]
        if m not in unk: d['Mass_flowrate']=TR[m]
        M.add_point(f,n,**d)
    pt('hi','Water',PB,'H_hi','mh'); pt('ho','Water',PB,'H_ho','mh'); pt('ci','R245fa',PEV,'H_ci','mc'); pt('co','R245fa',PEV,'H_co','mc')
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'hi','ho','ci','co',div_N=N); M.Component['E']=hx
    b=snap(M); t=time.perf_counter(); res='SOLVED'; msg=''
    try:
        with contextlib.redirect_stdout(io.StringIO()): hx.Cal()
    except Exception as e: res='refused'; msg=str(e)[:150].replace('\n',' ')
    a=snap(M)
    diff={k:{x:(b[k][x],a[k][x]) for x in b[k] if b[k][x]!=a[k].get(x)} for k in b}
    diff={k:v for k,v in diff.items() if v}
    print(label,sorted(unk),res,f'{time.perf_counter()-t:.1f}s',msg); print('   changed:',diff or 'none',flush=True)
# mid-solve failures (specification not consistent with the anchor)
trial({'H_hi','H_ho','mc'},label='mid-solve?')
trial({'H_hi','H_ci','mh'},label='mid-solve?',N=10)
trial({'H_hi','H_ho'},label='control 2-unknown')
