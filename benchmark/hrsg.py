import warnings
warnings.simplefilter('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
T0=300.0
def gt():
    # compressor
    h1=P('H','T',T0,'P',1e5,'Air'); s1=P('S','T',T0,'P',1e5,'Air')
    h2=h1+(P('H','P',12e5,'S',s1,'Air')-h1)/0.86
    h3=P('H','T',1400,'P',12e5,'Air'); s3=P('S','T',1400,'P',12e5,'Air')
    h4=h3-0.89*(h3-P('H','P',1e5,'S',s3,'Air'))
    T4=P('T','H',h4,'P',1e5,'Air')
    return dict(Wc=100*(h2-h1),Wt=100*(h3-h4),Q=100*(h3-h2),T4=T4)
G=gt()
def plant(p_bar,Tstack_C,ppt=15.0,N=200,mode='Evaporator'):
    Tin=G['T4']
    M=ThermodynamicModel()
    hs=P('H','P',0.008e6,'Q',0,'Water'); ss=P('S','P',0.008e6,'Q',0,'Water')
    hp=hs+(P('H','P',p_bar*1e5,'S',ss,'Water')-hs)/0.80
    M.add_point('Air','g_in',P=1e5,T=Tin,Mass_flowrate=100.0)
    M.add_point('Air','g_out',P=1e5,T=Tstack_C+273.15)
    M.add_point('Water','w_in',P=p_bar*1e5,H=hp)
    M.add_point('Water','w_out',P=p_bar*1e5,T=773.15)
    hx=HeatExchanger(M,'HRSG',ppt,'Evaporator',True,'g_in','g_out','w_in','w_out',div_N=N)
    hx.Cal()
    m=M.Point['w_in'].Mass_flowrate
    h5=P('H','T',773.15,'P',p_bar*1e5,'Water'); s5=P('S','T',773.15,'P',p_bar*1e5,'Water')
    h6s=P('H','P',0.008e6,'S',s5,'Water'); h6=h5-0.88*(h5-h6s)
    Wst=m*(h5-h6); Wp=m*(hp-hs)
    Wnet=G['Wt']-G['Wc']+Wst-Wp
    return dict(m=m,Q=hx.Q/1e6,mind=hx.min_dT,eta=100*Wnet/G['Q'],Wnet=Wnet/1e6,hx=hx,M=M,Wst=Wst/1e6,Tfw=M.Point['w_in'].T-273.15)
if __name__=='__main__':
    print({k:(v/1e6 if k!='T4' else v-273.15) for k,v in G.items()})
    for p,T in ((60,120),(40,150)):
        r=plant(p,T); print(p,T,{k:round(v,4) for k,v in r.items() if isinstance(v,float)})
