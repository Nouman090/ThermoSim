"""Condenser: R245fa (hot, condensing) vs water (cold). State built to be consistent with a cold-end anchor
T_ci = T_ho - 8 K, and with the minimum approach at that end (verified below)."""
import warnings, io, contextlib, signal, time, numpy as np
warnings.filterwarnings('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
PH=3e5; PC=2e5
from scipy.optimize import brentq
T_hi=343.15; T_ho=308.15; mh=2.0; T_ci=T_ho-8.0
H_hi=P('H','T',T_hi,'P',PH,'R245fa'); H_ho=P('H','T',T_ho,'P',PH,'R245fa'); H_ci=P('H','T',T_ci,'P',PC,'Water')
Q=mh*(H_hi-H_ho); H_dew=P('H','P',PH,'Q',1,'R245fa'); T_dew=P('T','P',PH,'Q',1,'R245fa')
Q_dew=mh*(H_dew-H_ho)                       # duty from the cold end up to the dew point
# choose mc so the approach at the dew point is ALSO 8 K  -> two pinch points (cold end + dew point)
mc=brentq(lambda m: T_dew-P('T','H',H_ci+Q_dew/m,'P',PC,'Water')-8.0, 5, 200)
H_co=H_ci+Q/mc; T_co=P('T','H',H_co,'P',PC,'Water')
f=np.linspace(0,1,4001); prof=[P('T','H',H_ho+x*(H_hi-H_ho),'P',PH,'R245fa')-P('T','H',H_ci+x*(H_co-H_ci),'P',PC,'Water') for x in f]
print(f'truth: Q={Q/1e3:.1f} kW mc={mc:.4f} kg/s T_co={T_co-273.15:.2f} C  dT cold end={prof[0]:.3f}  dT dew={T_dew-P("T","H",H_ci+Q_dew/mc,"P",PC,"Water"):.3f}  min over profile={min(prof):.3f} K')
TR=dict(H_hi=H_hi,H_ho=H_ho,H_ci=H_ci,H_co=H_co,mh=mh,mc=mc)
def _a(*x): raise TimeoutError
signal.signal(signal.SIGALRM,_a)
combos=[{'H_ci','H_hi','mh'},{'H_ci','H_hi','mc'},{'H_ci','H_co','mh'},{'H_ci','H_co','mc'},
        {'H_ho','H_hi','mh'},{'H_ho','H_hi','mc'},{'H_ho','H_co','mh'},{'H_ho','H_co','mc'}]
for unk in combos:
    M=ThermodynamicModel()
    def pt(n,f_,Pr,h,m):
        d=dict(P=Pr)
        if h not in unk: d['H']=TR[h]
        if m not in unk: d['Mass_flowrate']=TR[m]
        M.add_point(f_,n,**d)
    pt('hi','R245fa',PH,'H_hi','mh'); pt('ho','R245fa',PH,'H_ho','mh'); pt('ci','Water',PC,'H_ci','mc'); pt('co','Water',PC,'H_co','mc')
    hx=HeatExchanger(M,'C',8.0,'Condenser',None,'hi','ho','ci','co'); M.Component['C']=hx
    t=time.perf_counter(); signal.alarm(120)
    try:
        with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('always'); hx.Cal()
        signal.alarm(0)
        got=dict(H_hi=M.Point['hi'].H,H_ho=M.Point['ho'].H,H_ci=M.Point['ci'].H,H_co=M.Point['co'].H,mh=M.Point['hi'].Mass_flowrate,mc=M.Point['ci'].Mass_flowrate)
        e=max(abs(got[n]-TR[n])/abs(TR[n]) for n in unk)
        print(f"Condenser {{{','.join(sorted(unk))}}}: SOLVED {time.perf_counter()-t:4.1f}s maxrel={e:.1e} minDT={hx.min_dT:.3f} warns={len(w)}",flush=True)
    except TimeoutError: print(f"Condenser {{{','.join(sorted(unk))}}}: TIMEOUT",flush=True)
    except Exception as ex:
        signal.alarm(0); print(f"Condenser {{{','.join(sorted(unk))}}}: REFUSED {time.perf_counter()-t:4.1f}s {str(ex)[:200].replace(chr(10),' ')}",flush=True)
