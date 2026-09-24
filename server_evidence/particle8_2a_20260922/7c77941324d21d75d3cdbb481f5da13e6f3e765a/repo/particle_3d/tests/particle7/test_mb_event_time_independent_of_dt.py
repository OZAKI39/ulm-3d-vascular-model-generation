from particle_3d.injection_population import *
from particle_3d.particle7_cases import SONOVUE
import numpy as np
def run(dt):
 s=InjectionScheduler(PopulationSource(SONOVUE),LinearProfile([0],[1e-12])); rows=[]
 for t in np.arange(1,round(.25/dt)+1)*dt: rows.extend(s.through(float(t)))
 return rows
def test_exact():
 a,b,c=run(.25),run(.125),run(.0625)
 assert a==b==c and sum(e['species']=='MB' for e in a)==2
 assert len({e['scheduled_time_s'] for e in a})>2000
