"""Matched no-F2 single-bubble diagnostic, not an altered multibubble engine."""
import time
import numpy as np,h5py
from .common import *
from .single import integrate
from .static import force_context
from scripts.straight_tube_exp1.f10af.field import direct_force

def main():
 begin=time.perf_counter();verify();states=['near','mid','far','random_0','random_1','random_2'];write(BASE/'F1_SHADOW_PLAN.json',dict(states=states,duration_s=.15,cache='coarse',purpose='Same fixed initial IDs and field, independent F1-only dynamics versus full F1/F2/contact',identity=identity()));rows=[]
 with h5py.File(BASE/'INITIAL_STATES.h5') as h:
  for state in states:rows.append(integrate('F1_shadow_'+state,2500000,0,h['AG32_'+state][:],cache='coarse',duration=.15))
 source=[]
 for f in [2400000,2500000,2600000]:
  c,model,_,_,_=force_context(f)
  for j,name in enumerate(['Y-','Y+','X-','X+']):
   a=np.zeros(4);a[j]=1000;P,G,F,valid=direct_force(model,c['response'],np.zeros((1,3)),a,np.zeros(4));source.append(dict(frequency_Hz=f,source=name,pressure_Pa=P,gradient_Pa_m=G,force_N=F,response_valid=valid,scope='HOMOGENEOUS_WATER_PROXY_SPECIFIC'))
 write(BASE/'F1_SHADOW_RESULTS.json',dict(rows=rows,source_frequency_diagnostic=source,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
