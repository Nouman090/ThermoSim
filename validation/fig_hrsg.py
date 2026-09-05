import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, warnings
warnings.simplefilter('ignore')
from ThermoSim import *
from ThermoSim.plotting import CyclePlotter

def cc(P_hp, T_stack):
    M=ThermodynamicModel(); M.set_dead_state()
    rp,m_air,P_cond=12.0,100.0,0.008e6
    M.add_point('Air','g1',P=1e5,T=300,Mass_flowrate=m_air)
    M.add_point('Air','g2',P=1e5*rp)
    M.add_point('Air','g3',P=1e5*rp,T=1400,Mass_flowrate=m_air)
    M.add_point('Air','g4',P=1e5); M.add_point('Air','g5',P=1e5,T=T_stack)
    M.add_point('water','s1',P=P_hp,T=773.15); M.add_point('water','s2',P=P_cond)
    M.add_point('water','s3',P=P_cond,Q=0);    M.add_point('water','s4',P=P_hp)
    Compressor(M,'Compressor','g1','g2',n_isen=0.86,Calculate=True)
    Turbine(M,'GasTurbine','g3','g4',n_isen=0.89,Calculate=True)
    HeatExchanger(M,'Combustor',PPT=5,HEX_type='SimpleHEX',HeatAdded=True,
        Hot_In_state=None,Hot_Out_state=None,Cold_In_state='g2',Cold_Out_state='g3',Calculate=True)
    Pump(M,'Pump','s3','s4',n_isen=0.80,Calculate=True)
    HeatExchanger(M,'HRSG',PPT=15,HEX_type='double_pipe',HeatAdded=None,
        Hot_In_state='g4',Hot_Out_state='g5',Cold_In_state='s4',Cold_Out_state='s1',Calculate=True)
    Turbine(M,'SteamTurbine','s1','s2',n_isen=0.88,Calculate=True)
    HeatExchanger(M,'Condenser',PPT=5,HEX_type='SimpleHEX',HeatAdded=False,
        Hot_In_state='s2',Hot_Out_state='s3',Cold_In_state=None,Cold_Out_state=None,Calculate=True)
    M.Component['Pump'].Cal()
    M.ModelSummary(verbose=False)
    return M

cases = [(40e5, 423.15, 'feasible'), (60e5, 393.15, 'infeasible')]
fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0), sharey=False)

for ax, (P, Ts, tag) in zip(axes, cases):
    M = cc(P, Ts); p = CyclePlotter(M)
    p.plot_hex_profile('HRSG', div_N=200, save_csv=True, csv_path=f'/tmp/{tag}.csv')
    d = pd.read_csv(f'/tmp/{tag}.csv')
    Q = d['Heat_Transferred_kW']/1e3
    ax.plot(Q, d['Hot_Temperature_C'],  '-',  color='#c0392b', lw=2.0, label='exhaust gas')
    ax.plot(Q, d['Cold_Temperature_C'], '-',  color='#1f6fb4', lw=2.0, label='water / steam')
    mn = d['Delta_T_K'].min(); i = d['Delta_T_K'].idxmin()
    ax.plot([Q[i], Q[i]], [d['Cold_Temperature_C'][i], d['Hot_Temperature_C'][i]],
            color='#333', lw=1.1, ls='--')
    ax.annotate(f'min $\\Delta T$ = {mn:.1f} K', xy=(Q[i], (d['Cold_Temperature_C'][i]+d['Hot_Temperature_C'][i])/2),
                xytext=(10, -14), textcoords='offset points', fontsize=9,
                color='#c0392b' if mn < 0 else '#2e7d32', fontweight='bold', va='center')
    eta = M.Efficiency
    ax.set_title(f'({"a" if tag=="feasible" else "b"})  {P/1e5:.0f} bar, {Ts-273.15:.0f} °C stack\n'
                 f'$\\eta$ = {eta:.2f} %   —   {tag}', fontsize=10)
    ax.set_xlabel('cumulative duty [MW]'); ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8.5, loc='upper left')
axes[0].set_ylabel('temperature [°C]')
fig.tight_layout()
fig.savefig('fig_hrsg_pinch.png', dpi=300, bbox_inches='tight')
print('saved')
