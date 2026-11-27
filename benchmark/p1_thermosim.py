import time, warnings
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger

def build(pev_bar=15.0, ppt=8.0, heat_added=False, Tbrine=403.15):
    h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
    h4=h3+(P('H','P',pev_bar*1e5,'S',s3,'R245fa')-h3)/0.75
    M=ThermodynamicModel()
    M.add_point('Water','b_in', Mass_flowrate=10.0, T=Tbrine, P=5e5)
    M.add_point('Water','b_out', P=5e5)
    M.add_point('R245fa','w_in', P=pev_bar*1e5, H=h4)
    M.add_point('R245fa','w_out', P=pev_bar*1e5, Q=1)
    hx=HeatExchanger(M,'EVAP',ppt,'Evaporator',heat_added,'b_in','b_out','w_in','w_out')
    M.Component['EVAP']=hx
    return M,hx

for ha in (False,True,None):
    try:
        M,hx=build(heat_added=ha)
        t=time.perf_counter(); 
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            hx.Cal()
        dt=time.perf_counter()-t
        print(f'HeatAdded={ha}: OK  m_wf={M.Point["w_in"].Mass_flowrate:.4f} kg/s  '
              f'T_b_out={M.Point["b_out"].T-273.15:.2f} C  min_dT={hx.min_dT}  Q={hx.Q}  time={dt:.3f}s  warnings={[str(x.message)[:80] for x in w]}')
    except Exception as e:
        print(f'HeatAdded={ha}: FAIL {type(e).__name__}: {str(e)[:200]}')
