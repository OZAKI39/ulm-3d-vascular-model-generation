"""Bounded full-physics fallback, no LF-driven rejection or early stopping."""
import time
import numpy as np
from scipy.stats import qmc
from .common import *
from .run import task,batch

def main():
    begin=time.perf_counter();assert read(BASE/'RESTRICTED_LF_ADMISSION.json')['PASS']
    # Design fixed independently of responses. No claim of MF scientific admission.
    points=qmc.Sobol(3,scramble=True,seed=462060).random_base2(3);tasks=[]
    for i,u in enumerate(points):
        amp=500+1000*u[0];family='X_PAIR' if i<4 else 'XY_INPHASE';a=[0,0,amp,amp] if i<4 else [amp]*4;p=[0,0,0,float(2*np.pi*u[1])] if i<4 else [0,float(2*np.pi*u[2]),0,float(2*np.pi*u[1])]
        tasks.append(task(f'sobol_{i}',family,duration=.02,amplitudes=a,phases=p))
    write(BASE/'SOBOL_PLAN.json',dict(tasks=tasks,seed=462060,criterion='Structure factor gain and actual pattern gates, not latent surrogate score',physical_duration_s=.02,automatic_pruning=False,identity=identity()))
    rows=batch(tasks);write(BASE/'SOBOL_RESULTS.json',dict(rows=rows,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
