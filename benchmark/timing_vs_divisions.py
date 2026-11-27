import warnings, io, contextlib, time, statistics, json
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
import ThermoSim.model as tm
from ThermoSim import ThermodynamicModel, HeatExchanger
calls=[0]; orig=tm.ThermodynamicModel.Prop
def cnt(self,*a,**k): calls[0]+=1; return orig(self,*a,**k)
tm.ThermodynamicModel.Prop=cnt
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa'); h4=h3+(P('H','P',15e5,'S',s3,'R245fa')-h3)/0.75
def run(N,three):
    M=ThermodynamicModel(); M.add_point('Water','b_in',Mass_flowrate=10.0,T=403.15,P=5e5); M.add_point('Water','b_out',P=5e5)
    M.add_point('R245fa','w_in',P=15e5,H=h4); M.add_point('R245fa','w_out',**(dict(P=15e5) if three else dict(P=15e5,Q=1)))
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out',div_N=N); M.Component['E']=hx
    calls[0]=0; t=time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()): hx.Cal()
    return time.perf_counter()-t, calls[0], M.Point['w_in'].Mass_flowrate
out=[]
for three in (False,True):
    for N in (20,50,100,200):
        ts=[run(N,three) for _ in range(3)]
        out.append(dict(three=three,N=N,t=statistics.median(x[0] for x in ts),calls=ts[0][1],m=ts[0][2]))
        print(out[-1],flush=True)
json.dump(out,open('timing_vs_divisions.json','w'),indent=1)
