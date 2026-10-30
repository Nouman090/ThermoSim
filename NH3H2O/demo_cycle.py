"""Worked example: single-effect NH3-H2O absorption chiller with a rectifier.

Exercises the flash interface end to end -- (P,Q), (P,T), (P,H), (P,S) -- and
checks that mass, species and energy balances close. They close to machine
precision because the VLE is solved from the same Gibbs function that supplies
the enthalpies; pairing Ziegler-Trepp enthalpies with an external bubble/dew
correlation would leave a residual of several kW here.
"""
import warnings
from scipy.optimize import brentq
import nh3h2o as aw

warnings.filterwarnings('ignore')

P_hi, P_lo = 15.0, 4.0        # bar
T_gen, T_abs = 400.0, 310.0   # K
y_r = 0.995                   # NH3 mass fraction leaving the rectifier
m_s = 1.0                     # kg/s of strong solution

z_s = brentq(lambda z: aw.bubble_T(P_lo, z) - T_abs, 0.01, 0.99, xtol=1e-12)
st1 = aw.state(z=z_s, P=P_lo, Q=0)                 # absorber outlet
st2 = aw.state(z=z_s, P=P_hi, S=st1.s)             # isentropic pump
W_pump = m_s * (st2.h - st1.h)

gen = aw.state(z=z_s, P=P_hi, T=T_gen)             # generator flash
x_w = gen.x
st_weak = aw.state(z=x_w, P=P_hi, Q=0)
st_vap = aw.state(z=y_r, P=P_hi, Q=1)              # rectifier outlet

f = (z_s - x_w) / (y_r - x_w)                      # species balance
m_v, m_w = f * m_s, (1 - f) * m_s
Q_genset = m_v * st_vap.h + m_w * st_weak.h - m_s * st2.h

st_cond = aw.state(z=y_r, P=P_hi, Q=0)
Q_cond = m_v * (st_vap.h - st_cond.h)
st_ev_in = aw.state(z=y_r, P=P_lo, H=st_cond.h)    # isenthalpic throttle
st_ev_out = aw.state(z=y_r, P=P_lo, Q=1)
Q_evap = m_v * (st_ev_out.h - st_ev_in.h)

st_weak_lo = aw.state(z=x_w, P=P_lo, H=st_weak.h)
Q_abs = m_v * st_ev_out.h + m_w * st_weak_lo.h - m_s * st1.h
COP = Q_evap / (Q_genset + W_pump)

print('SINGLE-EFFECT NH3-H2O ABSORPTION CHILLER')
print(f'  P_hi = {P_hi:g} bar, P_lo = {P_lo:g} bar, T_gen = {T_gen:.0f} K, T_abs = {T_abs:.0f} K')
print(f'  strong z = {z_s:.5f}   weak x = {x_w:.5f}   refrigerant y = {y_r:.4f}')
print(f'  circulation ratio m_s/m_v = {1/f:.3f}')
print()
print(f'  {"state":<22s}{"T [K]":>9s}{"P [bar]":>9s}{"h [kJ/kg]":>11s}{"s [kJ/kg/K]":>13s}{"Q":>8s}')
for lbl, st in [('1 absorber out', st1), ('2 pump out', st2),
                ('  generator 2-phase', gen), ('3 rectified vapour', st_vap),
                ('4 condenser out', st_cond), ('5 evaporator in', st_ev_in),
                ('6 evaporator out', st_ev_out), ('7 weak, high P', st_weak),
                ('8 weak, low P', st_weak_lo)]:
    print(f'  {lbl:<22s}{st.T:9.2f}{st.P:9.3f}{st.h:11.2f}{st.s:13.5f}{st.Q:8.4f}')
print()
for lbl, v in [('Q_gen (gen+rect)', Q_genset), ('Q_cond', Q_cond),
               ('Q_evap', Q_evap), ('Q_abs', Q_abs), ('W_pump', W_pump)]:
    print(f'  {lbl:<18s}= {v:10.4f} kW')
print(f'  {"COP":<18s}= {COP:10.4f}')
print(f'  evaporator glide  = {st_ev_out.T - st_ev_in.T:9.2f} K '
      f'({st_ev_in.T:.2f} -> {st_ev_out.T:.2f} K)')

res = (Q_genset + W_pump + Q_evap) - (Q_cond + Q_abs)
scale = max(abs(Q_genset), abs(Q_evap))
print('\nBALANCE CHECKS')
print(f'  energy  in - out            : {res:12.4e} kW  ({100*abs(res)/scale:.2e} % of Q_gen)')
print(f'  total mass  m_v + m_w - m_s : {m_v + m_w - m_s:12.4e} kg/s')
print(f'  NH3 species                 : {m_v*y_r + m_w*x_w - m_s*z_s:12.4e} kg/s')
print(f'  throttle 4->5 isenthalpic   : {abs(st_ev_in.h - st_cond.h):12.4e} kJ/kg')
print(f'  throttle 7->8 isenthalpic   : {abs(st_weak_lo.h - st_weak.h):12.4e} kJ/kg')
print(f'  pump 1->2 isentropic        : {abs(st2.s - st1.s):12.4e} kJ/kg/K')
