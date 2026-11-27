import warnings, io, contextlib, json, itertools
warnings.simplefilter('ignore')
import CoolProp.CoolProp as CP
from ThermoSim import ThermodynamicModel, Turbine, Pump, HeatExchanger
from ThermoSim.plotting import CyclePlotter
import matplotlib; matplotlib.use('Agg')
def kH(P,k): 
    hf=CP.PropsSI('H','P',P,'Q',0,'R134a'); hg=CP.PropsSI('H','P',P,'Q',1,'R134a'); return hf+k*(hg-hf)
def run(Pev=16.0,Tgi=73.33,Tgo=54.44,Tcw=4.44,Te_in=12.85,kap=1.04,Pc=4.39,eta=0.80,Pw=2e5,N=300,sat_insert=True):
    M=ThermodynamicModel(); M.set_dead_state(T0=Tcw+273.15,P0=101325)
    M.add_point('R134a','1',P=Pev*1e5,H=kH(Pev*1e5,kap)); M.add_point('R134a','2',P=Pc*1e5)
    M.add_point('R134a','3',P=Pc*1e5,Q=0); M.add_point('R134a','4',P=Pev*1e5,T=Te_in+273.15)
    M.add_point('water','g1',P=Pw,T=Tgi+273.15,Mass_flowrate=33.39); M.add_point('water','g2',P=Pw,T=Tgo+273.15)
    M.add_point('water','c1',P=Pw,T=Tcw+273.15,Mass_flowrate=101.68); M.add_point('water','c2',P=Pw)
    with contextlib.redirect_stdout(io.StringIO()):
        HeatExchanger(M,'Evaporator',PPT=0.5,HEX_type='Evaporator',HeatAdded=True,Hot_In_state='g1',Hot_Out_state='g2',Cold_In_state='4',Cold_Out_state='1',Calculate=True)
        Turbine(M,'Turbine','1','2',n_isen=eta,Calculate=True)
        HeatExchanger(M,'Condenser',PPT=0.5,HEX_type='Condenser',HeatAdded=False,Hot_In_state='2',Hot_Out_state='3',Cold_In_state='c1',Cold_Out_state='c2',Calculate=True)
        Pump(M,'Pump','3','4',n_isen=0.75,Calculate=True)
        M.ModelSummary(verbose=False)
        p=CyclePlotter(M); de=p.plot_hex_profile('Evaporator',div_N=N); dc=p.plot_hex_profile('Condenser',div_N=N)
    return dict(m=M.Point['1'].Mass_flowrate,Tc2=M.Point['c2'].T-273.15,dEv=de,dCo=dc,Wt=M.Component['Turbine'].work/1e3,Qe=M.Component['Evaporator'].Q/1e3)
if __name__=='__main__':
    out={}
    for name,kw in dict(plant16=dict(),ips1695=dict(Pev=16.95),ips1695_g2=dict(Pev=16.95,Tgo=54.94)).items():
        for N in (50,100,300,500):
            r=run(N=N,**kw); out[f'{name}_{N}']=r; print(name,N,{k:round(v,4) for k,v in r.items()},flush=True)
    json.dump(out,open('chena_base.json','w'),indent=1)
