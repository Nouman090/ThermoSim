import time, warnings, io, contextlib
from CoolProp.CoolProp import PropsSI as P
from ThermoSim import ThermodynamicModel, HeatExchanger

PEV=15e5
h3=P('H','T',313.15,'Q',0,'R245fa'); s3=P('S','T',313.15,'Q',0,'R245fa')
h4=h3+(P('H','P',PEV,'S',s3,'R245fa')-h3)/0.75

def run(label, setup, strict=False):
    M=ThermodynamicModel()
    hx=setup(M, strict)
    buf=io.StringIO(); t=time.perf_counter()
    with warnings.catch_warnings(record=True) as w, contextlib.redirect_stdout(buf):
        warnings.simplefilter('always')
        try:
            hx.Cal(); status='returned'; err=''
        except Exception as e:
            status='EXCEPTION'; err=f'{type(e).__name__}: {e}'
    dt=time.perf_counter()-t
    print(f'--- {label} [{status}] {dt:.2f}s')
    if err: print('   exception :',err[:400])
    for x in w: print('   warning   :',str(x.message)[:400])
    out=buf.getvalue().strip()
    if out: print('   stdout    :',out[:400].replace('\n',' | '))
    try:
        print('   m_wf =',M.Point['w_in'].Mass_flowrate,' T_b_out =',(M.Point['b_out'].T-273.15) if M.Point['b_out'].T else None,
              ' T_w_out =',(M.Point['w_out'].T-273.15) if M.Point['w_out'].T else None,' min_dT =',hx.min_dT)
    except Exception as e: print('   (state read failed:',e,')')

def pts(M, bout, wout, bflow=10.0, wflow=None):
    M.add_point('Water','b_in', Mass_flowrate=bflow, T=403.15, P=5e5)
    M.add_point('Water','b_out', **bout)
    M.add_point('R245fa','w_in', Mass_flowrate=wflow, P=PEV, H=h4)
    M.add_point('R245fa','w_out', **wout)

# P2: all terminal states fixed -> pinch is only a CHECK. Brine out 90 C forces a crossing.
def p2(M,strict):
    pts(M, dict(P=5e5,T=363.15), dict(P=PEV,Q=1))
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out',strict_pinch=strict); M.Component['E']=hx; return hx
run('P2 crossing, warn mode', p2, strict=False)
run('P2 crossing, strict_pinch=True', p2, strict=True)

# P3: three unknowns {H_ho, H_co, m_c}; hot-end anchor + internal pinch (double pinch)
def p3(M,strict):
    pts(M, dict(P=5e5), dict(P=PEV))
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out',strict_pinch=strict); M.Component['E']=hx; return hx
run('P3 three unknowns {H_ho,H_co,m_c}', p3)

# P4a: case 15 {m_h, m_c}  (both flows unknown)
def p4a(M,strict):
    pts(M, dict(P=5e5,T=363.15), dict(P=PEV,Q=1), bflow=None)
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx; return hx
run('P4a unsupported: {m_h,m_c} unknown', p4a)

# P4b: case 1 {H_hi,H_ho}: hot side states unknown, cold side + both flows known
def p4b(M,strict):
    M.add_point('Water','b_in', Mass_flowrate=10.0, P=5e5)
    M.add_point('Water','b_out', P=5e5)
    M.add_point('R245fa','w_in', Mass_flowrate=4.7347, P=PEV, H=h4)
    M.add_point('R245fa','w_out', P=PEV, Q=1)
    hx=HeatExchanger(M,'E',8.0,'Evaporator',False,'b_in','b_out','w_in','w_out'); M.Component['E']=hx; return hx
run('P4b hot side unknown: {H_hi,H_ho}', p4b)
