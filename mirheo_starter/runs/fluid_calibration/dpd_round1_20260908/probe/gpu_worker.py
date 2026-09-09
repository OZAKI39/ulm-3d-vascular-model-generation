"""Actual small periodic liquid worker; invoked ONLY by the bounded MPI runner.

Adapted API wiring from Mirheo tests/flow/double_poiseuille.py and
tests/stress/pressure.py, Copyright 2019 ETH Zurich, MIT (MIRHEO_LICENSE.txt).
New initialization, safe chunk protocol and moment statistics are local code.
This file is self-contained so its exact snapshot can execute in a run folder.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import time
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec',required=True)
    args=parser.parse_args()  # --help ends before importing Mirheo/CUDA
    if os.environ.get('OMPI_COMM_WORLD_SIZE')!='2':
        raise RuntimeError('Worker requires authorized launcher with one compute + one postprocess rank')
    spec=json.loads(Path(args.spec).read_text());directory=Path.cwd()
    rank=int(os.environ['OMPI_COMM_WORLD_RANK']);is_compute=rank==0
    import mirheo as mir
    dt=spec['dt_star'];domain=tuple(spec['domain_star']);cand=spec['candidate'];task=spec['task']
    u=mir.Mirheo((1,1,1),domain,debug_level=3,log_filename='mirheo',no_splash=True)
    pv=mir.ParticleVectors.ParticleVector('pv',mass=spec['m_star'])
    u.registerParticleVector(pv,mir.InitialConditions.Uniform(number_density=task['n_star']))
    if is_compute:
        count=len(pv.getCoordinates())
        v=np.random.default_rng(spec['initial_velocity_seed']).normal(0,math.sqrt(spec['kBT_star']/spec['m_star']),(count,3))
        v-=v.mean(axis=0)
        pv.setVelocities(v.tolist())
    dpd=mir.Interactions.Pairwise('dpd',rc=spec['rc_star'],kind='DPD',a=cand['a'],gamma=cand['gamma'],kBT=spec['kBT_star'],power=cand['power'],stress=True,stress_period=dt/2)
    u.registerInteraction(dpd);u.setInteraction(dpd,pv,pv)
    vv=mir.Integrators.VelocityVerlet_withPeriodicForce('vv',force=task['force_star'],direction='x') if task['force_star'] else mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv);u.setIntegrator(vv,pv)
    every=spec['native_stats_every']
    u.registerPlugins(mir.Plugins.createStats('stats',every=every,filename='native_stats.csv'))
    u.registerPlugins(mir.Plugins.createVirialPressurePlugin('pressure',pv,lambda r:1.0,(1.,1.,1.),every,'pressure'))
    u.registerPlugins(mir.Plugins.createParticleChecker('checker',check_every=every))
    nb=int(round(domain[1]/spec['y_bin_star']));chunk=spec['snapshot_every'];done=0;index=0
    momfile=profilefile=None
    if is_compute:
        momfile=(directory/'moments.csv').open('x',buffering=1)
        profilefile=(directory/'profile_samples.csv').open('x',buffering=1)
        moments=csv.writer(momfile);profiles=csv.writer(profilefile)
        moments.writerow(['step','time_star','N','mass_star','kBT_raw_star','kBT_COM_star','kBT_thermal_star','mean_vx','mean_vy','mean_vz','max_speed_star','compute_wall_s','snapshot_wall_s'])
        profiles.writerow(['step','time_star','bin','y_star','count','ux','uy','uz','sum_v2'])
    started=time.monotonic();compute_total=0.;snapshot_total=0.
    while True:
        control=directory/'control'/f'chunk_{index:07d}.json'
        if is_compute:
            next_steps=0 if (directory/'STOP_REQUESTED').exists() else min(chunk,spec['steps']-done)
            tmp=control.with_suffix('.tmp')
            with tmp.open('x') as f:json.dump({'steps':next_steps},f)
            tmp.rename(control)
        else:
            # Both ranks read the same compute-authored stop decision, preventing
            # mismatched numbers of u.run() calls. Parent watchdog bounds waits.
            while not control.exists():time.sleep(.001)
        next_steps=json.loads(control.read_text())['steps']
        if not next_steps:break
        before=time.monotonic();u.run(next_steps,dt=dt);compute_total+=time.monotonic()-before;done+=next_steps
        if is_compute:
            before=time.monotonic()
            pos=np.asarray(pv.getCoordinates(),dtype=float);vel=np.asarray(pv.getVelocities(),dtype=float)
            if not np.isfinite(pos).all() or not np.isfinite(vel).all():raise RuntimeError('NONFINITE_PARTICLE_STATE')
            N=len(vel);mean=vel.mean(axis=0);raw=spec['m_star']*np.sum(vel**2)/(3*N)
            com=spec['m_star']*np.sum((vel-mean)**2)/(3*(N-1))
            ids=np.minimum(nb-1,np.floor(np.mod(pos[:,1],domain[1])/domain[1]*nb).astype(int))
            counts=np.bincount(ids,minlength=nb)
            sums=np.column_stack([np.bincount(ids,weights=vel[:,j],minlength=nb) for j in range(3)])
            av=np.divide(sums,counts[:,None],out=np.zeros_like(sums),where=counts[:,None]>0)
            v2=np.bincount(ids,weights=np.sum(vel**2,axis=1),minlength=nb)
            # Flow T removes instantaneous local mean in each y bin, with 3(N-B)
            # DOF. Residual within-bin shear biases upward; assess bin refinement
            # before production. Equilibrium uses only COM subtraction.
            thermal=spec['m_star']*np.sum((vel-av[ids])**2)/(3*(N-np.count_nonzero(counts))) if task['kind']=='flow' else com
            for j in range(nb):profiles.writerow([done,done*dt,j,(j+.5)*domain[1]/nb,int(counts[j]),*av[j],v2[j]])
            snap=time.monotonic()-before;snapshot_total+=snap
            moments.writerow([done,done*dt,N,N*spec['m_star'],raw,com,thermal,*mean,float(np.max(np.linalg.norm(vel,axis=1))),compute_total,snapshot_total])
        index+=1
    if is_compute:
        momfile.close();profilefile.close()
        with (directory/'worker_completion.json').open('x') as f:
            json.dump({'steps':done,'time_star':done*dt,'loop_wall_s':time.monotonic()-started,'compute_wall_s':compute_total,
                       'snapshot_wall_s':snapshot_total,'compute_s_per_step':compute_total/max(done,1),
                       'cooperative_stop':(directory/'STOP_REQUESTED').exists(),'initial_velocity_seed':spec['initial_velocity_seed'],
                       'dpd_native_pair_seed':42424242,'number_of_coordinators_per_rank':1},f,indent=2,allow_nan=False)
    del u


if __name__=='__main__':main()
