import warnings, io, contextlib
warnings.simplefilter('ignore')
import CoolProp.CoolProp as CP
from ThermoSim import ThermodynamicModel, HeatExchanger
def kH(P,k):
    hf=CP.PropsSI('H','P',P,'Q',0,'R134a'); hg=CP.PropsSI('H','P',P,'Q',1,'R134a'); return hf+k*(hg-hf)
for Pev,ppt in ((16.0,2.4017),(16.95,0.3068)):
    M=ThermodynamicModel(); P_=Pev*1e5
    M.add_point('R134a','1',P=P_,H=kH(P_,1.04)); M.add_point('R134a','4',P=P_,T=12.85+273.15)
    M.add_point('water','g1',P=2e5,T=73.33+273.15,Mass_flowrate=33.39); M.add_point('water','g2',P=2e5)
    hx=HeatExchanger(M,'E',ppt,'Evaporator',True,'g1','g2','4','1',div_N=200); M.Component['E']=hx
    with contextlib.redirect_stdout(io.StringIO()): hx.Cal()
    print(Pev,ppt,'brine exit',M.Point['g2'].T-273.15,'m',M.Point['1'].Mass_flowrate,'minDT',hx.min_dT)
