"""Re-test the 8 Evaporator 3-unknown combinations the anchor rules admit, using a state that IS consistent with
a hot-end anchor (P3 solution: T_hi=130, T_co=122.00, mc=4.1321, T_ho=106.19, 8 K at hot end and 8 K pinch)."""
import warnings, io, contextlib, itertools, signal, time, sys
warnings.filterwarnings('ignore')
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger
PEV=15e5; PB=5e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
H_ci=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75
H_hi=P('H','T',403.15,'P',PB,'Water'); H_co=P('H','T',395.15,'P',PEV,'R245fa'); mh=10.0; mc=4.1321
H_ho=H_hi-mc*(H_co-H_ci)/mh
TR=dict(H_hi=H_hi,H_ho=H_ho,H_ci=H_ci,H_co=H_co,mh=mh,mc=mc)
def _a(*x): raise TimeoutError
signal.signal(signal.SIGALRM,_a)
HEX=sys.argv[1] if len(sys.argv)>1 else 'Evaporator'
combos=[{'H_hi','H_ho','mh'},{'H_hi','H_ho','mc'},{'H_hi','H_ci','mh'},{'H_hi','H_ci','mc'},
        {'H_co','H_ho','mh'},{'H_co','H_ho','mc'},{'H_co','H_ci','mh'},{'H_co','H_ci','mc'}]
for unk in combos:
    M=ThermodynamicModel()
    def pt(n,f,Pr,h,m):
        d=dict(P=Pr)
        if h not in unk: d['H']=TR[h]
        if m not in unk: d['Mass_flowrate']=TR[m]
        M.add_point(f,n,**d)
    pt('hi','Water',PB,'H_hi','mh'); pt('ho','Water',PB,'H_ho','mh'); pt('ci','R245fa',PEV,'H_ci','mc'); pt('co','R245fa',PEV,'H_co','mc')
    hx=HeatExchanger(M,'E',8.0,HEX,False,'hi','ho','ci','co'); M.Component['E']=hx
    t=time.perf_counter(); signal.alarm(120)
    try:
        with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('always'); hx.Cal()
        signal.alarm(0)
        got=dict(H_hi=M.Point['hi'].H,H_ho=M.Point['ho'].H,H_ci=M.Point['ci'].H,H_co=M.Point['co'].H,mh=M.Point['hi'].Mass_flowrate,mc=M.Point['ci'].Mass_flowrate)
        e=max(abs(got[n]-TR[n])/abs(TR[n]) for n in unk)
        print(f"{HEX:11s} {{{','.join(sorted(unk))}}}: SOLVED {time.perf_counter()-t:4.1f}s maxrel={e:.1e} minDT={hx.min_dT:.3f} warns={len(w)}",flush=True)
    except TimeoutError: print(f"{HEX:11s} {{{','.join(sorted(unk))}}}: TIMEOUT",flush=True)
    except Exception as ex:
        signal.alarm(0); print(f"{HEX:11s} {{{','.join(sorted(unk))}}}: REFUSED {time.perf_counter()-t:4.1f}s {str(ex)[:210].replace(chr(10),' ')}",flush=True)
