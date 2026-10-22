"""Four formal controls, one cached acoustic basis, identical initial populations."""
import copy
import json
from pathlib import Path
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mevo_mb.coupling.experiment import run_particles, flow_from_config
from mevo_mb.io.artifacts import new_run as particle_run
from mevo_mb.utils.config import load_config
from scripts.stage3.common import new_run, atomic_json, finish, watermark


def main():
    c=load_config('configs/stage3/cylinder.yaml')
    history=json.loads(Path('results/stage3/history/STATE.json').read_text())
    assert history['status']=='PASS' and Path('docs/HISTORY_FORCE_SENSITIVITY.md').is_file()
    basis=Path(json.loads(Path('results/stage3/acoustic/STATE.json').read_text())['basis'])
    out=new_run('cylinder',dict(config=c,scenarios=['C1','C2','C3','C4'],target_transverse_radius_m=80e-6),['Stage3 scoped physics gates'])
    rows=[];positions={};summaries={};runs={}
    for name,primary,secondary,pattern in [('C1',False,False,'baseline'),('C2',True,False,'baseline'),('C3',True,True,'baseline'),('C4',True,True,'shifted')]:
        cc=copy.deepcopy(c);cc['particles'].update(primary_enabled=primary,secondary_enabled=secondary)
        if name=='C1':
            for p in cc['pzt']:p['enabled']=False
        run=particle_run(out,name,cc);runs[name]=str(run)
        try:
            summary=run_particles(cc,basis,run,pattern)
        except Exception as exc:
            atomic_json(out/'FAILURE.json',dict(scenario=name,error=repr(exc),run=str(run)))
            raise
        with h5py.File(run/'particles/trajectories.h5') as h:
            x=h['positions_m'][:];t=h['time_s'][:];flow=h['flow_velocity_m_s'][:];free=h['free_total_velocity_m_s'][:]
            primary_force=h['primary_force_N'][:];secondary_force=h['secondary_force_N'][:];local_pressure=abs(h['local_pressure_complex_Pa'][:])
        lateral=np.linalg.norm(x[-1,:,1:]-x[0,:,1:],axis=1)
        acoustic_speed=np.linalg.norm(free-flow,axis=-1);flow_speed=np.linalg.norm(flow,axis=-1)
        edges=np.linspace(-.0002,.0002,17)
        density0,_,_=np.histogram2d(x[0,:,1],x[0,:,2],bins=(edges,edges))
        density1,_,_=np.histogram2d(x[-1,:,1],x[-1,:,2],bins=(edges,edges))
        # Areal projected concentration; no inference of injected bulk concentration.
        cell_area=np.diff(edges)[:,None]*np.diff(edges)[None,:]
        np.savez_compressed(run/'cross_section_concentration.npz',edges_m=edges,initial_count=density0,final_count=density1,
            initial_areal_density_m2=density0/cell_area,final_areal_density_m2=density1/cell_area)
        row=dict(scenario=name,run=str(run),mean_lateral_displacement_m=float(lateral.mean()),
            lateral_displacement_percentiles_m={str(p):float(np.percentile(lateral,p)) for p in [10,50,90,95,100]},
            acoustic_drift_over_flow_mean=float(np.mean(acoustic_speed/np.maximum(flow_speed,1e-30))),
            acoustic_drift_over_flow_definition='Unconstrained force/drag drift magnitude over local flow magnitude, sampled at accepted states',
            initial_target_region_fraction=float(np.mean(np.linalg.norm(x[0,:,1:],axis=1)<80e-6)),
            final_target_region_fraction=float(np.mean(np.linalg.norm(x[-1,:,1:],axis=1)<80e-6)),
            target_definition='distance from cylinder axis < 80 um, fixed before runs',
            local_pressure_min_max_Pa=[float(local_pressure.min()),float(local_pressure.max())],
            mean_primary_force_magnitude_N=float(np.linalg.norm(primary_force,axis=-1).mean()),
            mean_secondary_force_magnitude_N=float(np.linalg.norm(secondary_force,axis=-1).mean()),
            initial=summary['initial'],final=summary['final'],pair_regime=summary['pair_regime'],
            NEAR_CONTACT_PHYSICS_APPROXIMATED=summary['NEAR_CONTACT_PHYSICS_APPROXIMATED'],timings=summary['timings'],
            validity_status=summary['validity_status'])
        if name=='C1':
            expected=x[0][None,:,:]+t[:,None,None]*flow_from_config(cc).velocity(x[0],0)[None,:,:]
            row['analytic_trajectory_error_m']=float(np.max(np.linalg.norm(x-expected,axis=-1)))
        rows.append(row);positions[name]=x;summaries[name]=summary
        atomic_json(out/'ROWS.json',rows)
    delta=np.linalg.norm(positions['C3'][-1]-positions['C4'][-1],axis=1)
    initial_equal=bool(np.array_equal(positions['C3'][0],positions['C4'][0]))
    ctrl3=json.loads((Path(runs['C3'])/'acoustic/CONTROL.json').read_text())
    ctrl4=json.loads((Path(runs['C4'])/'acoustic/CONTROL.json').read_text())
    only_phase=ctrl3['source_amplitudes_Pa']==ctrl4['source_amplitudes_Pa'] and initial_equal
    tests=dict(completed_four=all(q['status']=='COMPLETE' for q in summaries.values()),
        flow_only_analytic=rows[0]['analytic_trajectory_error_m']<1e-10,
        primary_lateral_migration=rows[1]['mean_lateral_displacement_m']>1e-8,
        secondary_changes_distribution=bool(np.linalg.norm(positions['C3'][-1]-positions['C2'][-1])>1e-8),
        phase_changes_distribution=bool(delta.mean()>1e-8),phase_only_controls=only_phase,
        no_penetration=all(q['final']['minimum_surface_gap_m']>=-c['particles']['contact_tolerance_m']*1.01 for q in rows))
    fig,axs=plt.subplots(1,4,figsize=(14,4),layout='constrained')
    for ax,(name,x) in zip(axs,positions.items()):
        ax.scatter(x[0,:,1]*1e6,x[0,:,2]*1e6,s=9,c='gray',alpha=.4,label='initial')
        ax.scatter(x[-1,:,1]*1e6,x[-1,:,2]*1e6,s=12,label='final')
        ax.add_patch(plt.Circle((0,0),80,fill=False,linestyle='--',color='green'))
        ax.set(title=name,xlabel='y (um)',ylabel='z (um)',aspect='equal');ax.legend()
    watermark(fig,'SCOPED MINIMAL MODEL; CONTACT IS GEOMETRIC; NO HISTORY FORCE')
    fig.savefig(out/'C1_C4_distribution.png',dpi=240);fig.savefig(out/'C1_C4_distribution.pdf');plt.close(fig)
    finish(out,dict(status='PASS' if all(tests.values()) else 'FAIL',tests=tests,rows=rows,
        C3_C4_mean_particle_position_difference_m=float(delta.mean()),
        C3_C4_cross_section_centroid_difference_m=float(np.linalg.norm(positions['C3'][-1,:,1:].mean(axis=0)-positions['C4'][-1,:,1:].mean(axis=0))),
        scope='Water cylinder, 64 monodisperse 1.78 um bubbles, 2.5 MHz, 2 ms. Engineering migration/aggregation trends; no validated transient displacement, near-field physics or experimental calibration.'))

if __name__=='__main__':main()
