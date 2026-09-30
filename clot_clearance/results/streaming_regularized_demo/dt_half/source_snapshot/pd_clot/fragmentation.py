"""Load-driven cyclic failure and continued PD fragment transport.

This entry point extends the project; it does not change the original runner.
Every point is retained, every force-bearing bond comes from damage.g, and the
only constrained positions are the original anchored particles.
"""
import argparse,csv,hashlib,json,time,traceback
from pathlib import Path
import numpy as np
from .geometry import make_cloud
from .mechanics import timestep_limit
from .damage import particle_damage
from .fragment_mechanics import prepare_supported_shape,evaluate_supported
from .fragment_topology import BondDamage,FragmentTracker,surface_particles
from .fragment_fluid import FragmentFluid,relaxation_half_step
from .output import write_json


def run(c, output, quiet=False):
    out=Path(output);out.mkdir(parents=True,exist_ok=False);write_json(out/'CONFIG.json',c)
    cloud=make_cloud(c['clot']);x=cloud.X.copy();v=np.zeros_like(x)
    dt=c['simulation']['dt_s'];limit=timestep_limit(cloud,c)
    if not 0<dt<=limit:raise ValueError(f'Unsafe explicit timestep {dt}; estimate {limit}')
    dc=c['damage'];s=c['simulation'];period=1/s['representative_frequency_Hz']
    steps_per_cycle=int(round(period/dt))
    if abs(steps_per_cycle*dt-period)>1e-12:raise ValueError('dt must divide the representative cycle')
    if c['material']['model']!='neo_hookean':raise ValueError('Unsupported material')
    damage=BondDamage(len(cloud.pairs),dc);fluid=FragmentFluid(c)
    tracker=FragmentTracker(cloud,c['transport']['x_clearance_m'])
    mass=cloud.volume*c['clot']['density_kg_m3']
    shape=prepare_supported_shape(cloud,damage.g,c['safety'])
    attached,components,stats=tracker.update(damage.g,x,v,0)
    exposed,bulk_ratio=surface_particles(cloud,damage.g,c['surface'])
    initial_exposed=exposed.copy();states=[];history=[];component_history=[];N=0;t=0.;tick=0
    impulse=np.zeros_like(x);cum_impulse=np.zeros_like(x)
    events={name:None for name in ['first_bond_failure','first_detachment','first_fragment','first_clearance_crossing','first_reduced_rank']}
    root=Path(__file__).resolve().parents[1]
    source_files=[*root.glob('pd_clot/*.py'),*root.glob('vendor/*.py')]
    write_json(out/'IDENTITY.json',dict(model='existing_NOSB_PD_failure_transport_extension',
        no_prescribed_cleavage=True,initial_positions='unaltered cloud.X',initial_velocity='zero',
        config_sha256=hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest(),
        source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        particle_count=len(x),bond_count=len(damage.g),dt_estimate_s=limit))
    vtk=out/'vtk';vtk.mkdir();start=time.perf_counter();minimum_J=1.;max_acceleration=0.

    def forces():
        mech=evaluate_supported(cloud,x,damage.g,c['material'],c['safety'],shape)
        surface_force,normals,traction=fluid.surface_force(cloud,x,damage.g,exposed,attached,t)
        acc=(mech[0]+surface_force)/mass[:,None];acc[cloud.fixed]=0
        if not np.isfinite(acc).all():raise FloatingPointError('Nonfinite acceleration')
        return mech,acc,normals,traction

    def record(step, recent, amplitude):
        import pyvista as pv
        mech,_,normals,traction=forces();D,Dmax=particle_damage(cloud,damage.g)
        row=dict(macro_step=step,represented_cycles=N,mechanical_time_s=t,
            mean_damage=float(D.mean()),maximum_damage=float(D.max()),
            active_bond_count=int(damage.active.sum()),broken_bond_count=int((~damage.active).sum()),
            broken_bond_fraction=float(np.mean(~damage.active)),particle_count=len(x),total_volume_m3=float(cloud.volume.sum()),
            maximum_displacement_m=float(np.linalg.norm(x-cloud.X,axis=1).max()),
            rank_counts=[int(np.sum(shape[3]==r)) for r in range(4)],
            surface_particle_count=int(exposed.sum()),new_surface_particle_count=int(np.sum(exposed&~initial_exposed)),
            maximum_cycle_stretch_amplitude=float(amplitude.max()),
            hydrodynamic_impulse_norm_kg_m_s=float(np.linalg.norm(impulse)),
            base_max_error_m=float(np.abs(x[cloud.fixed]-cloud.X[cloud.fixed]).max()),**stats)
        history.append(row);component_history.append(dict(cycles=N,mechanical_time_s=t,components=components))
        data=dict(x=x.copy(),v=v.copy(),integrity=damage.g.copy(),bond_D=damage.D.copy(),bond_active=damage.active.copy(),
            damage=D.copy(),fragment_id=tracker.labels.copy(),attached=attached.copy(),surface=exposed.copy(),
            rank=shape[3].copy(),cleared=tracker.cleared.copy(),hydro_impulse=impulse.copy(),cycles=N,
            recent_broken=recent.copy(),amplitude=amplitude.copy())
        states.append(data)
        poly=pv.PolyData(x);poly['reference_position_m']=cloud.X;poly['displacement']=x-cloud.X;poly['velocity']=v
        poly['damage']=D;poly['particle_damage']=D;poly['max_bond_damage']=Dmax;poly['fragment_id']=tracker.labels
        poly['attached_to_base']=attached.astype(np.uint8);poly['is_surface_particle']=exposed.astype(np.uint8)
        poly['FIXED_BASE']=cloud.fixed.astype(np.uint8);poly['cleared_material']=tracker.cleared.astype(np.uint8)
        poly['number_of_active_bonds']=np.bincount(cloud.pairs[damage.active].ravel(),minlength=len(x))
        poly['deformation_rank']=shape[3];poly['J_or_intrinsic_stretch']=np.nan_to_num(mech[2],nan=0.)
        poly['strain_energy_density']=mech[3];poly['stabilization_energy_density']=mech[4]
        poly['particle_volume_m3']=cloud.volume;poly['fluid_velocity_m_s']=fluid.velocity(x,t)
        poly['macro_hydrodynamic_impulse_kg_m_s']=impulse;poly['surface_normal']=normals
        poly['surface_traction_Pa']=traction
        poly['normal_traction_Pa']=np.sum(traction*normals,axis=1)[:,None]*normals
        poly['tangential_traction_Pa']=traction-poly['normal_traction_Pa']
        poly.save(vtk/f'particles_{step:04d}.vtp')
        bonds=pv.PolyData(x,lines=np.column_stack((np.full(len(cloud.pairs),2),cloud.pairs)).ravel())
        bonds.cell_data['integrity']=damage.g;bonds.cell_data['accumulated_D']=damage.D
        bonds.cell_data['active']=damage.active.astype(np.uint8);bonds.cell_data['recently_broken']=recent.astype(np.uint8)
        bonds.cell_data['measured_cycle_amplitude']=amplitude;bonds.save(vtk/f'bonds_{step:04d}.vtp')
        write_json(out/'PROGRESS.json',row)
        if not quiet:print(json.dumps({k:row[k] for k in ['macro_step','represented_cycles','mean_damage','broken_bond_fraction','detached_volume_fraction','detached_components','cleared_volume_fraction','maximum_fragment_travel_m','rank_counts']}),flush=True)

    try:
        record(0,np.zeros(len(damage.g),bool),np.zeros(len(damage.g)))
        with (out/'failure_ledger.jsonl').open('w') as ledger:
            for macro in range(1,s['number_of_macro_steps']+1):
                lo=np.full(len(damage.g),np.inf);hi=np.full(len(damage.g),-np.inf);impulse[:]=0
                warmup=s.get('warmup_cycles',1) if macro==1 else 0
                nsteps=steps_per_cycle*(s['representative_cycles']+warmup)
                mech,acc,_,_=forces()
                free=(~attached)&(~cloud.fixed)
                damping=np.where(attached,s['damping_s_inv'],s['detached_damping_s_inv'])
                decay=np.exp(-damping[:,None]*dt/2)
                for k in range(nsteps):
                    v*=decay
                    if c['transport']['enabled']:
                        updated=relaxation_half_step(v,fluid.velocity(x,t),free,dt/2,c['transport']['tau_h_s'])
                        impulse+=mass[:,None]*(updated-v);v=updated
                    v+=.5*dt*acc;x+=dt*v;x[cloud.fixed]=cloud.X[cloud.fixed];v[cloud.fixed]=0
                    tick+=1;t=tick*dt
                    if np.linalg.norm(x-cloud.X,axis=1).max()>c['safety']['maximum_displacement_m']:
                        raise FloatingPointError('Configured displacement guard exceeded')
                    radius=np.linalg.norm(x[:,1:3],axis=1)
                    if radius.max()>c['pipe']['radius_m']*c['safety']['maximum_radius_ratio']:
                        raise FloatingPointError('Particle left prescribed pipe cross-section; no contact model active')
                    mech,acc,_,_=forces()
                    peak=float(np.linalg.norm(acc,axis=1).max());max_acceleration=max(max_acceleration,peak)
                    if peak>c['safety']['maximum_acceleration_m_s2']:raise FloatingPointError('Acceleration guard exceeded')
                    if np.any(shape[3]>0):minimum_J=min(minimum_J,float(mech[2][shape[3]>0].min()))
                    v+=.5*dt*acc
                    if c['transport']['enabled']:
                        updated=relaxation_half_step(v,fluid.velocity(x,t),free,dt/2,c['transport']['tau_h_s'])
                        impulse+=mass[:,None]*(updated-v);v=updated
                    v*=decay;v[cloud.fixed]=0
                    if k>=warmup*steps_per_cycle:
                        stretch=np.linalg.norm(x[cloud.pairs[:,1]]-x[cloud.pairs[:,0]],axis=1)/cloud.length
                        lo=np.minimum(lo,stretch);hi=np.maximum(hi,stretch)
                amplitude=.5*(hi-lo);cum_impulse+=impulse
                ids,previous_D,increment=damage.advance(amplitude);N+=dc['DeltaN']
                recent=np.zeros(len(damage.g),bool);recent[ids]=True
                for bid in ids:
                    ledger.write(json.dumps(dict(cycles=N,bond_id=int(bid),particle_ids=cloud.pairs[bid].tolist(),
                        Q=float(amplitude[bid]),D_before=float(previous_D[bid]),D_after=float(damage.D[bid]),
                        increment=float(increment[bid]),D_break=dc.get('D_break',1.)))+'\n')
                ledger.flush()
                shape=prepare_supported_shape(cloud,damage.g,c['safety'])
                attached,components,stats=tracker.update(damage.g,x,v,N)
                exposed,bulk_ratio=surface_particles(cloud,damage.g,c['surface'])
                conditions=dict(first_bond_failure=bool(np.any(~damage.active)),first_detachment=bool(np.any(~attached)),
                    first_fragment=stats['total_components']>1,first_clearance_crossing=bool(tracker.cleared.any()),
                    first_reduced_rank=bool(np.any(shape[3]<3)))
                for name,condition in conditions.items():
                    if condition and events[name] is None:events[name]=dict(cycles=N,macro_step=macro,mechanical_time_s=t)
                record(macro,recent,amplitude)
        summary=dict(status='COMPLETED_SYNTHETIC_FRAGMENTATION_RUN',label=c['verification_label'],
            elapsed_wall_s=time.perf_counter()-start,events=events,final=history[-1],minimum_J_or_intrinsic_stretch=minimum_J,
            maximum_internal_and_surface_acceleration_m_s2=max_acceleration,
            all_particles_retained=len(x)==len(cloud.X),cleavage_commands_used=False,
            physical_predictivity_claim=False,maximum_fragment_travel_m=max(h['maximum_fragment_travel_m'] for h in history))
        write_json(out/'SUMMARY.json',summary);write_json(out/'history.json',history)
        write_json(out/'components.json',component_history);write_json(out/'lineage.json',tracker.lineage)
        write_json(out/'clearance_ledger.json',tracker.crossings)
        with (out/'history.csv').open('w',newline='') as f:
            keys=[k for k,val in history[0].items() if not isinstance(val,list)]
            w=csv.DictWriter(f,keys,extrasaction='ignore');w.writeheader();w.writerows(history)
        np.savez_compressed(out/'states.npz',X=cloud.X,pairs=cloud.pairs,volume=cloud.volume,fixed=cloud.fixed,
            **{key:np.array([state[key] for state in states]) for key in states[0]})
        for kind in ['particles','bonds']:
            entries='\n'.join(f'<DataSet timestep="{row["represented_cycles"]}" file="vtk/{kind}_{i:04d}.vtp"/>' for i,row in enumerate(history))
            (out/f'{kind}.pvd').write_text('<?xml version="1.0"?><VTKFile type="Collection" version="0.1"><Collection>'+entries+'</Collection></VTKFile>')
        return summary
    except Exception as error:
        write_json(out/'FAILED.json',dict(error=repr(error),cycles=N,mechanical_time_s=t,events=events,traceback=traceback.format_exc()))
        np.savez_compressed(out/'failed_state.npz',X=cloud.X,x=x,v=v,pairs=cloud.pairs,integrity=damage.g,bond_D=damage.D)
        write_json(out/'history_partial.json',history)
        raise


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(json.loads(args.config.read_text()),args.output)


if __name__=='__main__':main()
