"""Slow-time erosion with fixed, independently configured mechanical windows.

Derived integration logic from the frozen regularization runner. Constitutive
mechanics, surface force and energy ledger are reused without editing that path.
"""
import argparse,csv,hashlib,json,shutil,time,traceback,warnings
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from ..geometry import make_cloud
from ..mechanics import timestep_limit
from ..damage import particle_damage
from ..fragment_mechanics import prepare_supported_shape
from ..fragment_topology import surface_particles
from ..fragment_fluid import relaxation_half_step
from ..output import write_json
from .fast_mechanics import evaluate_supported
from .flow import ConsistentFragmentFluid,ManufacturedStreamingProvider
from .resolved_flow import ResolvedStreamingFieldProvider
from .energy import kinetic,potential
from .erosion_energy import ErosionEnergyLedger
from .diagnostics import classify,localization
from .erosion import StreamingErosionDamage,tensile_strain,SettlingMonitor

ROOT=Path(__file__).resolve().parents[2]


def validate(c):
    if c['damage']['mode']!='STREAMING_EROSION' or not c['damage']['enabled']:
        raise ValueError('Independent STREAMING_EROSION mode required')
    if c['damage'].get('D_break')!=1.:raise ValueError('D_break must equal 1')
    if c['streaming']['oscillatory_amplitude']!=0. or c['streaming']['mean_fraction']!=1.:
        raise ValueError('STREAMING_EROSION v1 requires mean-only streaming')
    if c['streaming']['provider'] not in ['manufactured','resolved_field']:raise ValueError('Unsupported mean field')
    if c['transport']['mode']!='relaxation':raise ValueError('This erosion step reuses relaxation transport only')
    ec=c['streaming_erosion'];mc=c['erosion_mechanics'];dt=c['simulation']['dt_s']
    StreamingErosionDamage(1,ec)
    if not np.isfinite(ec['erosion_exposure_step_s']) or ec['erosion_exposure_step_s']<=0:raise ValueError('Positive erosion step required')
    steps=ec['number_of_erosion_steps']
    if isinstance(steps,bool) or int(steps)!=steps or steps<1:raise ValueError('Positive integer erosion step count required')
    if not np.isfinite(dt) or dt<=0 or not np.isfinite(mc['relaxation_time_s']) or mc['relaxation_time_s']<=0:raise ValueError('Positive mechanical interval required')
    n=round(mc['relaxation_time_s']/dt)
    if n<2 or not np.isclose(n*dt,mc['relaxation_time_s'],rtol=0,atol=1e-14):raise ValueError('Mechanical dt must divide relaxation_time_s')
    if not 0<mc['settling_check_fraction']<=1:raise ValueError('Invalid settling fraction')
    if not np.isfinite(mc['maximum_relative_driver_drift']) or mc['maximum_relative_driver_drift']<0:raise ValueError('Invalid drift tolerance')
    return int(n)


def graph_state(cloud,g,x,v,ranks,c):
    """Exact graph at this state; IDs are state-local, not cycle-based lineage."""
    pairs=cloud.pairs[g>0]
    graph=coo_matrix((np.ones(2*len(pairs)),(pairs.ravel(),pairs[:,::-1].ravel())),shape=(len(x),len(x)))
    count,labels=connected_components(graph,directed=False)
    attached=np.isin(labels,np.unique(labels[cloud.fixed]));components=[]
    for label in range(count):
        ids=np.flatnonzero(labels==label);weights=cloud.volume[ids]
        components.append(dict(fragment_id=label,particle_ids=ids.tolist(),particle_count=len(ids),attached_to_base=bool(attached[ids[0]]),
                               volume_m3=float(weights.sum()),center_of_mass_m=np.average(x[ids],axis=0,weights=weights).tolist(),
                               velocity_of_center_of_mass_m_s=np.average(v[ids],axis=0,weights=weights).tolist()))
    category,components,stats=classify(cloud,g,ranks,attached,labels,components,c['fragment_resolution'])
    total=float(cloud.volume.sum());stats.update(attached_volume_fraction=float(cloud.volume[attached].sum()/total),
        detached_volume_fraction=float(cloud.volume[~attached].sum()/total),component_count=count)
    return attached,labels,category,components,stats


def aggregate_bond_values(cloud,values):
    weights=cloud.weight*np.asarray(values)
    top=np.bincount(cloud.pairs.ravel(),weights=np.repeat(weights,2),minlength=len(cloud.X))
    bottom=np.bincount(cloud.pairs.ravel(),weights=np.repeat(cloud.weight,2),minlength=len(cloud.X))
    return np.divide(top,bottom,out=np.zeros_like(top),where=bottom>0)


def run(c,output,quiet=False,provider=None):
    nsteps=validate(c);out=Path(output);out.mkdir(parents=True,exist_ok=False)
    write_json(out/'CONFIG.json',c)
    cloud=make_cloud(c['clot']);cloud.fixed=cloud.X[:,2]<c['clot']['origin_m'][2]+c['clot']['fixed_base_thickness_m']-cloud.spacing*1e-10
    if not cloud.fixed.any() or cloud.fixed.all():raise ValueError('Invalid discretized fixed slab')
    dt=c['simulation']['dt_s'];limit=timestep_limit(cloud,c)
    if dt>limit:raise ValueError(f'Unsafe explicit timestep {dt}, estimated bound {limit}')
    ec,mc=c['streaming_erosion'],c['erosion_mechanics'];de=ec['erosion_exposure_step_s']
    damage=StreamingErosionDamage(len(cloud.pairs),ec);x=cloud.X.copy();v=np.zeros_like(x)
    if provider is None:
        provider=(ManufacturedStreamingProvider(c) if c['streaming']['provider']=='manufactured' else
                  ResolvedStreamingFieldProvider(c['streaming']['field_path'],c['streaming']['metadata_path'],c['pipe']['viscosity_Pa_s'],cloud.X[cloud.exposed],out/'IMPORT_DIAGNOSTIC.json'))
    if not isinstance(provider,(ManufacturedStreamingProvider,ResolvedStreamingFieldProvider)):
        raise ValueError('v1 requires a verified manufactured mean field or static resolved provider')
    for field_time in [0.,mc['relaxation_time_s']]:
        sample=provider.sample(cloud.X,field_time)
        if sample.du_dt is not None and np.any(sample.du_dt!=0):raise ValueError('Non-static flow is not supported by this mean-only mode')
    fluid=ConsistentFragmentFluid(c,provider);mass=cloud.volume*c['clot']['density_kg_m3']
    shape=prepare_supported_shape(cloud,damage.g,c['safety'])
    attached,labels,category,components,graph_stats=graph_state(cloud,damage.g,x,v,shape[3],c)
    exposed,_=surface_particles(cloud,damage.g,c['surface']);tick=0;mechanical_time=0.;exposure=0.;macro=0
    history=[];states=[];settling_windows=[];component_history=[];failures=[];events={'first_bond_failure':None,'first_detachment':None}
    source=[*ROOT.glob('pd_clot/*.py'),*ROOT.glob('pd_clot/regularization/*.py'),*ROOT.glob('vendor/*.py')]
    for p in source:
        dest=out/'source_snapshot'/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    identity=dict(mode='STREAMING_EROSION',config_sha256=hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest(),
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source},particle_count=len(x),bond_count=len(cloud.pairs),
        mean_only=True,independent_exposure_clock=True,fixed_mechanical_steps_per_macro=nsteps,mechanical_dt_estimate_s=limit,
        driver='Endpoint tensile strain; no temporal stretch range',parameters='DEVELOPMENT / UNCALIBRATED PLACEHOLDERS',
        coupon_calibration_used=False,flow_provider=type(provider).__name__,component_ids='State-local exact connected-component labels',
        transport_clock='Mechanical integration windows only; not treatment-time transport')
    write_json(out/'IDENTITY.json',identity);(out/'vtk').mkdir()
    start=time.perf_counter();minimum_J=1.;maximum_acceleration=0.

    def strain():return tensile_strain(np.linalg.norm(x[cloud.pairs[:,1]]-x[cloud.pairs[:,0]],axis=1)/cloud.length)

    def forces():
        mech=evaluate_supported(cloud,x,damage.g,c['material'],c['safety'],shape)
        force,normals,traction=fluid.surface_force(cloud,x,damage.g,exposed,attached,mechanical_time)
        acc=(mech[0]+force)/mass[:,None];acc[cloud.fixed]=0.
        if not np.isfinite(acc).all():raise FloatingPointError('Nonfinite acceleration')
        return mech,acc,force,normals,traction

    ledger=ErosionEnergyLedger(sum(potential(cloud,forces()[0])));bond_energy_drop=np.zeros(len(damage.D))

    def damping(decay):
        before=kinetic(mass,v);v[:]*=decay;ledger.damping_dissipation_J+=before-kinetic(mass,v)

    def transport(free):
        if c['transport']['enabled'] and np.any(free):
            before=kinetic(mass,v)
            v[:]=relaxation_half_step(v,fluid.velocity(x,mechanical_time),free,dt/2,c['transport']['tau_h_s'])
            ledger.transport_work_J+=kinetic(mass,v)-before

    def record(settling,driver,previous_damage):
        import pyvista as pv
        mech,_,_,normals,traction=forces();pd,pmax=particle_damage(cloud,damage.g)
        elastic,stab=potential(cloud,mech);er=ledger.row(elastic,stab,kinetic(mass,v))
        if abs(er['relative_numerical_energy_residual'])>c['energy_audit']['maximum_relative_residual'] and abs(er['numerical_energy_residual_J'])>c['energy_audit']['minimum_absolute_tolerance_J']:
            raise FloatingPointError(f'Unexplained energy residual {er}')
        current=strain();active_eps=current[damage.active];particle_driver=aggregate_bond_values(cloud,driver)
        row=dict(macro_step=macro,mechanical_time_s=mechanical_time,erosion_exposure_time_s=exposure,
            mechanical_relaxation_time_s=mc['relaxation_time_s'],erosion_exposure_step_s=de,
            mean_damage=float(pd.mean()),maximum_damage=float(pd.max()),mean_bond_damage=float(damage.D.mean()),maximum_bond_damage=float(damage.D.max()),
            maximum_epsilon_plus=float(active_eps.max()) if len(active_eps) else 0.,mean_epsilon_plus=float(active_eps.mean()) if len(active_eps) else 0.,
            maximum_used_erosion_rate_s_inv=float(damage.last_rate_s_inv.max()),broken_bond_count=int((~damage.active).sum()),
            active_bond_count=int(damage.active.sum()),particle_count=len(x),total_volume_m3=float(cloud.volume.sum()),
            maximum_displacement_m=float(np.linalg.norm(x-cloud.X,axis=1).max()),base_max_error_m=float(np.abs(x[cloud.fixed]-cloud.X[cloud.fixed]).max()),
            rank_counts=[int(np.sum(shape[3]==k)) for k in range(4)],**graph_stats,**er,
            **localization(cloud,pd,damage.g,attached,np.array(c['streaming']['bubble_center_m']),c['localization']['influence_radius_m']),
            **(settling or dict(settling_status='INITIAL_STATE_NOT_TESTED',driver_relative_drift=None)))
        history.append(row);component_history.append(dict(erosion_exposure_time_s=exposure,components=components))
        states.append(dict(x=x.copy(),v=v.copy(),integrity=damage.g.copy(),bond_D=damage.D.copy(),bond_active=damage.active.copy(),
            damage=pd,rank=shape[3].copy(),attached=attached.copy(),fragment_id=labels.copy(),resolution_class=category.copy(),
            epsilon_plus=current,erosion_driver_epsilon_plus=driver.copy(),particle_erosion_driver=particle_driver,
            erosion_rate_s_inv=damage.last_rate_s_inv.copy(),D_before_update=previous_damage.copy(),bond_stored_energy_drop_J=bond_energy_drop.copy(),
            mechanical_time_s=mechanical_time,erosion_exposure_time_s=exposure))
        poly=pv.PolyData(x)
        for key,value in dict(reference_position_m=cloud.X,displacement_m=x-cloud.X,velocity_m_s=v,damage=pd,maximum_bond_damage=pmax,
            particle_erosion_driver=particle_driver,fixed_base=cloud.fixed.astype(np.uint8),rank=shape[3],resolution_class=category,
            surface_traction_Pa=traction,surface_normal=normals,particle_volume_m3=cloud.volume).items():poly[key]=value
        poly.save(out/f'vtk/particles_{macro:04d}.vtp')
        bonds=pv.PolyData(x,lines=np.column_stack((np.full(len(cloud.pairs),2),cloud.pairs)).ravel())
        for key,value in dict(integrity=damage.g,D=damage.D,active=damage.active.astype(np.uint8),epsilon_plus=current,
                              used_erosion_driver=driver,rate_s_inv=damage.last_rate_s_inv).items():bonds.cell_data[key]=value
        bonds.save(out/f'vtk/bonds_{macro:04d}.vtp');write_json(out/'PROGRESS.json',row)
        if not quiet:print(json.dumps({k:row[k] for k in ['macro_step','mechanical_time_s','erosion_exposure_time_s','mean_damage','maximum_damage','broken_bond_count','detached_total_volume_m3','settling_status','driver_relative_drift']}),flush=True)

    try:
        record(None,np.zeros(len(damage.D)),damage.D.copy())
        first_check=max(0,int(np.ceil(nsteps*(1-mc['settling_check_fraction']))) - 1)
        for macro in range(1,int(ec['number_of_erosion_steps'])+1):
            monitor=SettlingMonitor(damage.active,mc['driver_absolute_floor']);window_start=mechanical_time
            free=(~attached)&(~cloud.fixed)
            decay=np.exp(-np.where(attached,c['simulation']['damping_s_inv'],c['simulation']['detached_damping_s_inv'])[:,None]*dt/2)
            mech,acc,fs,_,_=forces()
            for k in range(nsteps):
                damping(decay);transport(free);v[:]+=.5*dt*acc;dx=dt*v;x[:]+=dx
                x[cloud.fixed]=cloud.X[cloud.fixed];v[cloud.fixed]=0.;tick+=1;mechanical_time=tick*dt
                if np.linalg.norm(x-cloud.X,axis=1).max()>c['safety']['maximum_displacement_m']:raise FloatingPointError('Displacement guard exceeded')
                if np.linalg.norm(x[:,1:3],axis=1).max()>c['pipe']['radius_m']*c['safety']['maximum_radius_ratio']:raise FloatingPointError('Particle left pipe cross-section')
                old_fs=fs;mech,acc,fs,_,_=forces();ledger.external_work_J+=float(np.sum(.5*(old_fs+fs)*dx))
                maximum_acceleration=max(maximum_acceleration,float(np.linalg.norm(acc,axis=1).max()))
                if maximum_acceleration>c['safety']['maximum_acceleration_m_s2']:raise FloatingPointError('Acceleration guard exceeded')
                if np.any(shape[3]>0):minimum_J=min(minimum_J,float(mech[2][shape[3]>0].min()))
                v[:]+=.5*dt*acc;transport(free);damping(decay);v[cloud.fixed]=0.
                if k>=first_check:monitor.add(mechanical_time,mechanical_time-window_start,strain(),kinetic(mass,v))
            settling=monitor.summary(mc['maximum_relative_driver_drift'])
            settling_windows.append(dict(macro_step=macro,erosion_exposure_start_s=(macro-1)*de,erosion_exposure_end_s=macro*de,
                                         mechanical_window_start_s=window_start,summary=settling,samples=monitor.rows))
            if settling['settling_status']=='MECHANICAL_STEADY_STATE_NOT_ESTABLISHED':
                warnings.warn('MECHANICAL_STEADY_STATE_NOT_ESTABLISHED: retaining endpoint strain with explicit diagnostic flag',RuntimeWarning)
            driver=strain();before=sum(potential(cloud,mech));ids,old,increment=damage.advance(driver,de);exposure=macro*de
            shape=prepare_supported_shape(cloud,damage.g,c['safety']);new_mech=evaluate_supported(cloud,x,damage.g,c['material'],c['safety'],shape)
            if np.any(shape[3]>0):minimum_J=min(minimum_J,float(new_mech[2][shape[3]>0].min()))
            # Only a diagnostic allocation by damage increment, not a fracture calibration.
            increases_before=ledger.damage_energy_increase_count
            bond_energy_drop+=ledger.damage_drop(before,sum(potential(cloud,new_mech)),increment)
            if ledger.damage_energy_increase_count>increases_before:
                warnings.warn('ENERGY_DISSIPATION_NOT_ESTABLISHED: damage increased stored energy; signed change retained',RuntimeWarning)
            for bid in ids:failures.append(dict(erosion_exposure_time_s=exposure,mechanical_time_s=mechanical_time,bond_id=int(bid),
                particle_ids=cloud.pairs[bid].tolist(),D_before=float(old[bid]),D_after=float(damage.D[bid]),
                epsilon_plus=float(driver[bid]),erosion_rate_s_inv=float(damage.last_rate_s_inv[bid])))
            attached,labels,category,components,graph_stats=graph_state(cloud,damage.g,x,v,shape[3],c)
            exposed,_=surface_particles(cloud,damage.g,c['surface'])
            for name,condition in [('first_bond_failure',bool(len(ids))),('first_detachment',bool(np.any(~attached)))]:
                if condition and events[name] is None:events[name]=dict(erosion_exposure_time_s=exposure,mechanical_time_s=mechanical_time,macro_step=macro)
            record(settling,driver,old)
        statuses=[h['settling_status'] for h in history[1:]]
        summary=dict(status='COMPLETED_STREAMING_EROSION_DEVELOPMENT_RUN',scientific_status='SOFTWARE / NUMERICAL DEVELOPMENT',
            parameters='DEVELOPMENT / UNCALIBRATED PLACEHOLDERS',physical_calibration='NOT PERFORMED',mesh_objectivity='NOT ASSESSED IN THIS STEP',
            wall_attachment_mesh_consistency='STILL OPEN',combined_erosion_fatigue='NOT IMPLEMENTED',real_microbubble_CFD='NOT USED IN THIS VERIFICATION',
            coupon_calibration_used=False,D_break=1.,final=history[-1],events=events,
            failure_observation='NO FAILURE / NO DETACHMENT OBSERVED' if events['first_bond_failure'] is None and events['first_detachment'] is None else 'EVENTS_OBSERVED',
            mechanical_settling_status='DRIVER_SETTLING_CHECK_PASSED' if all(s=='DRIVER_SETTLING_CHECK_PASSED' for s in statuses) else 'MECHANICAL_STEADY_STATE_NOT_ESTABLISHED',
            maximum_driver_relative_drift=max(h['driver_relative_drift'] for h in history[1:]),
            minimum_J_or_intrinsic_stretch=minimum_J,maximum_acceleration_m_s2=maximum_acceleration,elapsed_wall_s=time.perf_counter()-start,
            all_particles_retained=True,final_configuration_convention='Relaxation endpoint then damage update at fixed x; no extra post-update settling window',
            energy_scope='Actual mechanical path only; never multiplied by slow exposure',
            damage_energy_allocation='SIGNED stored-energy drop split proportional to D increment; not positive dissipation and not Gc',
            damage_energy_status=history[-1]['damage_energy_status'])
        for name,value in [('SUMMARY.json',summary),('history.json',history),('components.json',component_history),('settling_windows.json',settling_windows)]:write_json(out/name,value)
        np.savez_compressed(out/'states.npz',X=cloud.X,pairs=cloud.pairs,volume=cloud.volume,fixed=cloud.fixed,
                            **{key:np.array([s[key] for s in states]) for key in states[0]})
        with (out/'history.csv').open('w',newline='') as f:
            keys=list(dict.fromkeys(k for row in history for k,v in row.items() if not isinstance(v,list)))
            writer=csv.DictWriter(f,keys,extrasaction='ignore');writer.writeheader();writer.writerows(history)
        (out/'failure_ledger.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in failures))
        for kind in ['particles','bonds']:
            entries='\n'.join(f'<DataSet timestep="{h["erosion_exposure_time_s"]}" file="vtk/{kind}_{i:04d}.vtp"/>' for i,h in enumerate(history))
            (out/f'{kind}.pvd').write_text('<?xml version="1.0"?><VTKFile type="Collection" version="0.1"><Collection>'+entries+'</Collection></VTKFile>')
        return summary
    except Exception as exc:
        write_json(out/'FAILED.json',dict(error=repr(exc),traceback=traceback.format_exc(),macro_step=macro,mechanical_time_s=mechanical_time,erosion_exposure_time_s=exposure))
        np.savez_compressed(out/'failed_state.npz',X=cloud.X,x=x,v=v,bond_D=damage.D,integrity=damage.g)
        write_json(out/'history_partial.json',history);write_json(out/'settling_partial.json',settling_windows);raise


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(json.loads(a.config.read_text()),a.output)


if __name__=='__main__':main()
