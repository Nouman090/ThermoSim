import warnings; warnings.simplefilter('ignore')
import CoolProp.CoolProp as CP
# direct property-model sensitivity at the 16.95 bar pinch (bubble point), independent of ThermoSim
for P in (16.0e5,16.95e5):
    print(P/1e5,'Tsat HEOS',CP.PropsSI('T','P',P,'Q',0,'R134a')-273.15)
    for f in ('PR::R134a','SRK::R134a'):
        try: print('  ',f,CP.PropsSI('T','P',P,'Q',0,f)-273.15)
        except Exception as e: print('  ',f,'err',str(e)[:50])
for T in (54.44,60.0,73.33):
    a=CP.PropsSI('H','P',2e5,'T',T+273.15,'water'); b=CP.PropsSI('H','P',2e5,'T',T+273.15,'IF97::water')
    print('water h HEOS vs IF97 at',T,(b-a)/a*100,'%', 'cp diff', CP.PropsSI('C','P',2e5,'T',T+273.15,'IF97::water')-CP.PropsSI('C','P',2e5,'T',T+273.15,'water'))
