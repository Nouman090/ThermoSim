import warnings, sys
from ThermoSim import ThermodynamicModel, Turbine, Pump, HeatExchanger
def run(pev=15.0, ppt=8.0, n_t=0.80, n_p=0.75, verbose=False):
    M=ThermodynamicModel(); M.set_dead_state()
    M.add_point('R245fa','1',P=pev*1e5,Q=1); M.add_point('R245fa','2',P=2.506e5); M.add_point('R245fa','3',P=2.506e5,Q=0); M.add_point('R245fa','4',P=pev*1e5)
    M.add_point('water','b1',P=5e5,T=403.15,Mass_flowrate=10); M.add_point('water','b2',P=5e5)
    Pump(M,'Pump','3','4',n_isen=n_p,Calculate=True)
    HeatExchanger(M,'Evaporator',PPT=ppt,HEX_type='Evaporator',HeatAdded=True,Hot_In_state='b1',Hot_Out_state='b2',Cold_In_state='4',Cold_Out_state='1',Calculate=True)
    Turbine(M,'Turbine','1','2',n_isen=n_t,Calculate=True)
    HeatExchanger(M,'Condenser',PPT=5,HEX_type='SimpleHEX',HeatAdded=False,Hot_In_state='2',Hot_Out_state='3',Cold_In_state=None,Cold_Out_state=None,Calculate=True)
    M.Component['Pump'].Cal()
    M.ModelSummary(verbose=verbose)
    return M
if __name__=='__main__':
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always'); M=run(verbose=True)
    print([str(x.message)[:100] for x in w])
    print({k:v for k,v in vars(M).items() if not k.startswith('_') and isinstance(v,(int,float,str))})
    for c in M.Component.values(): print(c.ID, c.Ex_D, getattr(c,'W',None), getattr(c,'Q',None))
    for name,p in M.Point.items(): print(name, p.T-273.15, p.P/1e5, p.Mass_flowrate, getattr(p,'Ex',None))
