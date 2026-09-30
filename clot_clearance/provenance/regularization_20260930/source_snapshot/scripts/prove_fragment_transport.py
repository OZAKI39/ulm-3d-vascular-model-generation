"""Controlled replay from an ACTUAL detached state; no cleavage or resets.

Freeze its surviving bonds for a short force-attribution comparison. Both
branches start with identical saved x/v. The only branch difference is whether
the configured hydrodynamic relaxation force is enabled.
"""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pd_clot.geometry import make_cloud
from pd_clot.fragment_mechanics import prepare_supported_shape,evaluate_supported
from pd_clot.fragment_fluid import FragmentFluid,relaxation_half_step


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();run=args.run
    z=np.load(run/'states.npz');c=json.loads((run/'CONFIG.json').read_text());h=json.loads((run/'history.json').read_text())
    k=next(i for i in range(len(h)) if h[i]['detached_volume_fraction']>0)
    free=~z['attached'][k];g=make_cloud(c['clot']);b=z['integrity'][k]
    shape=prepare_supported_shape(g,b,c['safety']);provider=FragmentFluid(c)
    dt=c['simulation']['dt_s'];duration=1/c['simulation']['representative_frequency_Hz'];nstep=round(duration/dt)
    outcomes={}
    for enabled in [False,True]:
        x=z['x'][k].copy();v=z['v'][k].copy();initial=x.copy();t=h[k]['mechanical_time_s']
        impulse=np.zeros_like(x);mass=g.volume*c['clot']['density_kg_m3']
        decay=np.exp(-np.where(free,c['simulation']['detached_damping_s_inv'],c['simulation']['damping_s_inv'])[:,None]*dt/2)
        for step in range(nstep):
            acc=evaluate_supported(g,x,b,c['material'],c['safety'],shape)[0]/mass[:,None];acc[g.fixed]=0
            v*=decay
            if enabled:
                new=relaxation_half_step(v,provider.velocity(x,t),free,dt/2,c['transport']['tau_h_s']);impulse+=mass[:,None]*(new-v);v=new
            v+=.5*dt*acc;x+=dt*v;x[g.fixed]=g.X[g.fixed];v[g.fixed]=0;t+=dt
            acc=evaluate_supported(g,x,b,c['material'],c['safety'],shape)[0]/mass[:,None];acc[g.fixed]=0
            v+=.5*dt*acc
            if enabled:
                new=relaxation_half_step(v,provider.velocity(x,t),free,dt/2,c['transport']['tau_h_s']);impulse+=mass[:,None]*(new-v);v=new
            v*=decay;v[g.fixed]=0
        outcomes['hydro_on' if enabled else 'hydro_off']=dict(
            free_COM_displacement_m=np.average(x[free]-initial[free],axis=0,weights=g.volume[free]).tolist(),
            free_COM_velocity_m_s=np.average(v[free],axis=0,weights=g.volume[free]).tolist(),
            free_hydrodynamic_impulse_kg_m_s=impulse[free].sum(axis=0).tolist(),free_final_positions_m=x[free].tolist())
    difference=float(np.linalg.norm(np.array(outcomes['hydro_on']['free_COM_displacement_m'])-outcomes['hydro_off']['free_COM_displacement_m']))
    result=dict(status='PASS' if difference>1e-8 else 'FAIL',source_cycles=int(z['cycles'][k]),free_particle_ids=np.flatnonzero(free).tolist(),
        replay_duration_s=duration,initial_positions_and_velocities_identical=True,bonds_frozen_for_control_only=True,
        no_kinematic_separation=True,no_particle_velocity_reset=True,position_difference_due_to_hydrodynamic_force_m=difference,**outcomes)
    (run/'TRANSPORT_FORCE_ATTRIBUTION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if result['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
