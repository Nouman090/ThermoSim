import warnings, json; warnings.simplefilter('ignore')
import io, contextlib, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from orc import run
def Wn(**kw):
    with contextlib.redirect_stdout(io.StringIO()):
        M=run(**kw)
    t=M.Component['Turbine']; p=M.Component['Pump']
    return (abs(t.work)-abs(p.work))/1e3 if hasattr(t,'work') else None
W8=Wn(ppt=8.0); W4=Wn(ppt=4.0)
print('W(8K)',W8,'W(4K)',W4,'dW',W4-W8)
Wt80=Wn(n_t=0.80); Wt82=Wn(n_t=0.82); print('eta .80',Wt80,'.82',Wt82,'dW',Wt82-Wt80)
out=dict(W8=W8,W4=W4,dW_ppt=W4-W8,Wt80=Wt80,Wt82=Wt82,dW_eta=Wt82-Wt80)
# arc (midpoint) elasticity over the finite change, as the tool's formula
def arc(x0,x1,y0,y1): return (y1-y0)/(x1-x0)*((x0+x1)/2)/((y0+y1)/2)
out['arc_ppt']=arc(8,4,W8,W4); out['arc_eta']=arc(.80,.82,Wt80,Wt82)
# point elasticities, central differences
for h in (0.5,0.25,0.1):
    a,b=Wn(ppt=8+h),Wn(ppt=8-h); out[f'pt_ppt_h{h}']=(a-b)/(2*h)*8/W8
    out[f'grad_ppt_h{h}']=(a-b)/(2*h)
for h in (0.005,0.0025):
    a,b=Wn(n_t=.80+h),Wn(n_t=.80-h); out[f'pt_eta_h{h}']=(a-b)/(2*h)*.80/Wt80
# local at midpoint 6 K for comparison
a,b=Wn(ppt=6.25),Wn(ppt=5.75); W6=Wn(ppt=6.0); out['pt_ppt_6K']=(a-b)/0.5*6/W6
# log-log (constant-elasticity) secant
import math
out['loglog_ppt']=math.log(W4/W8)/math.log(4/8); out['loglog_eta']=math.log(Wt82/Wt80)/math.log(.82/.80)
for k,v in out.items(): print(k,round(v,5))
json.dump(out,open('elas.json','w'),indent=1)
