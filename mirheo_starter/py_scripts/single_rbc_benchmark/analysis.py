"""CPU analysis of stored evidence; never launches a solver."""
import math
import time
from pathlib import Path
import numpy as np
from .physics import read_off, geometry, unwrap, paths, groups, speedup
from .quality import mesh_intersections,leakage_probes,fluid_summary,local_flow_summary,feedback_comparison,dynamics_window,strict_comparison
from py_scripts.fluid_physics.common import read_json, write_json, atomic_state, sha256_file, fingerprint, now

def read_csv(path):
    if Path(path).stat().st_size==0:return np.array([],dtype=[])
    return np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True,dtype=None,encoding='utf8'))

def end_to_end_cost(execution_s,completion,analysis_s,case_creation_s=0.):
    # wall_preparation_s is contained in setup_s and must not be added twice.
    fields=('setup_s','relaxation_s','coupled_s','output_s','state_handoff_s')
    parts={k:float(completion.get(k,0.)) for k in fields}
    if min([execution_s,analysis_s,case_creation_s,*parts.values()])<0:raise ValueError('NEGATIVE_COST')
    if sum(parts.values())>execution_s+.01:raise ValueError('OVERLAPPING_PHASE_TIMES')
    return dict(process_s=execution_s,case_creation_s=case_creation_s,analysis_s=analysis_s,
                total_s=execution_s+case_creation_s+analysis_s,
                unassigned_process_overhead_s=max(0.,execution_s-sum(parts.values())),**parts)

def affine_response(path,area):
    a=read_csv(path);rows=[]
    for mode in np.unique(a['mode']):
        origin=a[(a['mode']==mode)&(a['epsilon']==0)];q=np.column_stack([origin[k] for k in ('x','y','z')]);baseforce=np.column_stack([origin[k] for k in ('fx','fy','fz')])
        for eps in np.unique(a['epsilon']):
            r=a[(a['mode']==mode)&(a['epsilon']==eps)];f=np.column_stack([r[k] for k in ('fx','fy','fz')])-baseforce
            derivative=np.column_stack([q[:,2],q[:,0]*0,q[:,0]*0]) if mode=='shear' else np.column_stack([q[:,0],q[:,0]*0,-q[:,2]/(1+eps)**2])
            force=-float(np.sum(f*derivative))
            rows.append(dict(mode=str(mode),epsilon=float(eps),generalized_restoring_force=force,apparent_modulus=force/(eps*area) if eps else None,force_l2=float(np.linalg.norm(f))))
    shear=[r['apparent_modulus'] for r in rows if r['mode']=='shear' and abs(abs(r['epsilon'])-.01)<1e-8]
    return dict(Gs_affine_effective=float(np.mean(shear)),definition='native total incremental force work under x += epsilon*z, divided by A0*epsilon; geometry-dependent effective response, not exact continuum shear modulus',rows=rows)

def periodic_poiseuille_basis(y,L,force):
    # u = force/nu * piecewise quadratic. Force is acceleration when m=1.
    y=np.asarray(y);q=np.where(y<L/2,y-L/4,y-3*L/4);sign=np.where(y<L/2,-1.,1.)
    return sign*force*((L/4)**2-q*q)/2

def viscosity(path,c):
    a=read_csv(Path(path)/'profiles.csv');a=a[(a['phase']=='shear')&(a['time_star']>=c['material_check']['drive_time']*(1-c['material_check']['fit_final_fraction']))]
    if len(a)==0:raise ValueError('NO_STATIONARY_MATERIAL_SAMPLES')
    times=np.unique(a['time_star']);estimates=[];residuals=[]
    for t in times:
        r=a[a['time_star']==t];b=periodic_poiseuille_basis(r['z'],c['material_check']['domain'][1],c['material_check']['force_per_particle']/c['dpd']['mass'])
        slope=float(np.dot(b,r['ux'])/np.dot(b,b));estimates.append(1/slope if slope>0 else np.nan);residuals.append(float(np.linalg.norm(r['ux']-b*slope)/max(np.linalg.norm(b*slope),1e-12)))
    # Average velocities before inversion; preserve per-frame scatter and block estimates.
    mean=np.array([a[a['bin']==b]['ux'].mean() for b in sorted(np.unique(a['bin']))]);y=np.array([a[a['bin']==b]['z'][0] for b in sorted(np.unique(a['bin']))]);b=periodic_poiseuille_basis(y,c['material_check']['domain'][1],c['material_check']['force_per_particle'])
    invnu=float(np.dot(b,mean)/np.dot(b,b));nu=1/invnu
    if not math.isfinite(nu) or nu<=0:raise ValueError('INVALID_MEASURED_VISCOSITY')
    return dict(nu=nu,mean_profile_relative_l2=float(np.linalg.norm(mean-b*invnu)/np.linalg.norm(b*invnu)),sample_count=len(times),frame_nu_std=float(np.nanstd(estimates)),
                method='known acceleration double-periodic Poiseuille; fit final half of drive; no Couette-viscosity inference',frame_nu=[float(x) if math.isfinite(x) else None for x in estimates])

def _analyze_record(directory,c):
    directory=Path(directory);execution=read_json(directory/'execution.json');complete=read_json(directory/'completion.json') if (directory/'completion.json').exists() else {}
    spec=read_json(directory/'spec.json') if (directory/'spec.json').exists() else {}
    result=dict(task=directory.name,directory=str(directory),execution=execution,completion=complete,frames=[],screen='NOT_EVALUATED',analysis_s=0.,invalid_saved_frames=[])
    t=time.perf_counter();mesh=directory/'reference.off'
    if not mesh.exists():mesh=paths(c)[0]/'common_reference.off'
    raw=directory/'vertices.csv'
    if raw.exists() and raw.stat().st_size>120 and mesh.exists():
        _,f=read_off(mesh);a=read_csv(raw)
        # Two native coordinators meet at the same global step, with distinct phases.
        # Preserve both snapshots instead of combining their vertex IDs.
        keys=list(dict.fromkeys(zip(a['phase'].tolist(),a['step'].tolist())))
        for phase,step in keys:
            r=a[(a['phase']==phase)&(a['step']==step)];r=r[np.argsort(r['vertex'])]
            try:
                if not np.array_equal(r['vertex'],np.arange(int(f.max())+1)):raise ValueError('BAD_SAVED_VERTEX_IDS')
                v=unwrap(np.column_stack([r[k] for k in ('x','y','z')]),f,[c['geometry']['periodic_length']]*2+[c['geometry']['gap']])
                if result['frames']:
                    previous=np.array(result['frames'][-1]['metrics']['center']);L=c['geometry']['periodic_length'];v[:,:2]+=L*np.rint((previous[:2]-v.mean(axis=0)[:2])/L)
                g=geometry(v,f,c['screen']['angle_degenerate_axis_ratio']);intersections=mesh_intersections(v,f)
                forces=np.column_stack([r[k] for k in ('fx','fy','fz')]);g['saved_force_rms']=float(np.sqrt(np.mean(np.sum(forces**2,axis=1))));g['saved_force_max']=float(np.linalg.norm(forces,axis=1).max());g['marker0']=v[0].tolist()
                result['frames'].append(dict(step=int(step),time_star=float(r['time_star'][0]),strain=float(r['strain'][0]),phase=str(r['phase'][0]),metrics=g,vertices=v.tolist(),intersection_check=intersections))
            except ValueError as error:result['invalid_saved_frames'].append(dict(phase=str(phase),step=int(step),error=str(error)))
        result['triangles']=f.tolist()
        first=result['frames'][0]['metrics'];metrics=[r['metrics'] for r in result['frames']]
        checks=dict(closed_single_cell=not result['invalid_saved_frames'] and all(g['closed'] and g['euler_characteristic']==2 for g in metrics),orientation_preserved=all(g['signed_volume']*first['signed_volume']>0 for g in metrics),
                    area_drift=max(abs(g['area']/first['area']-1) for g in metrics),volume_drift=max(abs(g['volume']/first['volume']-1) for g in metrics),
                    no_wall_crossing=all(g['bounds'][0][2]>0 and g['bounds'][1][2]<c['geometry']['gap'] for g in metrics),
                    deformation_observed=max(abs(g['D']-first['D']) for g in metrics)>1e-5,
                    strict_impermeability='NOT_VERIFIED',self_intersection_test='PASS_SAVED_NONADJACENT_TRIANGLES' if all(x['intersection_check']['nonadjacent_intersections']==0 for x in result['frames']) else 'FAILED',space_convergence='NOT_VERIFIED',
                    nonadjacent_intersections=sum(x['intersection_check']['nonadjacent_intersections'] for x in result['frames']))
        checks['geometry_screen']='PASS' if checks['closed_single_cell'] and checks['orientation_preserved'] and checks['no_wall_crossing'] and checks['nonadjacent_intersections']==0 and checks['area_drift']<=c['screen']['area_relative_drift'] and checks['volume_drift']<=c['screen']['volume_relative_drift'] else 'FAILED'
        checks['maximum_saved_vertex_force']=max(g['saved_force_max'] for g in metrics)
        checks['unsaved_force_limit_events']='NOT_INSTRUMENTED; HemoCell native 50 pN cap retained; sample maxima cannot establish that no intermediate limit event occurred.'
        result['checks']=checks
        result['leakage_probes']=leakage_probes(directory,result['frames'],f,c,{'prep_steps':spec.get('prep_steps',0),**complete})
        relaxed=[x for x in result['frames'] if x['phase']=='relaxation'];started=[x for x in result['frames'] if x['phase']=='shear' and abs(x['strain'])<1e-9]
        if relaxed and started:
            delta=np.array(relaxed[-1]['vertices'])-np.array(started[0]['vertices'])
            result['handoff_check']=dict(maximum_vertex_position_change=float(np.abs(delta).max()),source='Actual saved before/after meshes',velocities='Native same-process copy in worker; no independent velocity snapshot comparison saved')
    profile=directory/'profiles.csv'
    if profile.exists() and profile.stat().st_size>110:
        a=read_csv(profile);last=a[a['step']==max(a['step'])]
        result['last_profile']={k:last[k].tolist() for k in ('z','ux','density')}
    if complete.get('role',spec.get('role'))!='material':result['fluid']=fluid_summary(directory,c,read_csv)
    result['local_flow']=local_flow_summary(directory,read_csv)
    phases=complete
    if not complete and (directory/'timings.csv').exists():
        a=read_csv(directory/'timings.csv');phases=dict(relaxation_s=float(a[a['phase']=='relaxation']['compute_s'].sum()),coupled_s=float(a[a['phase']=='shear']['compute_s'].sum()),output_s=float(a['output_s'].sum())) if len(a) else {}
        result['partial_timing_note']='Completed timing chunks only; unassigned process time also includes the fatal/interrupted chunk and unrecorded setup. Actual native total step count after the last saved frame is unknown.'
    result['observed_last_saved_strain']=max([f['strain'] for f in result['frames']],default=None)
    result['observed_last_saved_step']=max([f['step'] for f in result['frames']],default=None)
    resource=directory/'host_resources.json';result['host_resources']=read_json(resource) if resource.exists() else None
    result['analysis_s']=time.perf_counter()-t
    case=directory/'case_creation.json'
    result['case_creation_measurement']='MEASURED' if case.exists() else 'NOT_RECORDED_FOR_EARLY_PREFLIGHT'
    result['costs']=end_to_end_cost(execution['elapsed_monotonic_s'],phases,result['analysis_s'],read_json(case)['wall_s'] if case.exists() else 0.)
    return result

def run_summary(directory,c):
    directory=Path(directory);base,_,_=paths(c)
    sealfile=directory/'output_sha256.json'
    key=fingerprint(dict(evidence=sha256_file(sealfile),config=c['_config_sha256'],
                         analysis=sha256_file(Path(__file__)),quality=sha256_file(Path(__file__).with_name('quality.py')),
                         geometry=sha256_file(Path(__file__).with_name('physics.py'))))
    cache=base/'analysis_cache'/directory.name/(key+'.json')
    if cache.exists():
        r=read_json(cache);r['analysis_cache_reused']=True;return r
    r=_analyze_record(directory,c);r['analysis_cache_reused']=False;r['analysis_cache_key']=key
    cache.parent.mkdir(parents=True,exist_ok=True);write_json(cache,r)
    with (base/'analysis_events.jsonl').open('a') as f:
        f.write(__import__('json').dumps(dict(task=directory.name,recorded_at=now(),analysis_s=r['analysis_s'],key=key))+'\n')
    return r


def collect(c):
    from py_scripts.solver_benchmark.workflow import verify_run
    base,runs,_=paths(c);all_runs=[]
    for p in sorted(runs.rglob('execution.json')):
        verify_run(p.parent)
        if p.parent.name.startswith(('compile','cmake')):
            all_runs.append(dict(task=p.parent.name,directory=str(p.parent),execution=read_json(p),completion={},frames=[],analysis_s=0.))
        else:all_runs.append(run_summary(p.parent,c))
    lookup={r['task']:r for r in all_runs};expected=[f'main_{engine}_{i}' for engine in ('hemocell','mirheo') for i in (1,2)]
    formal=[lookup[n] for n in expected if n in lookup]
    def completed(r):
        cp=r['completion'];coupling=cp.get('membrane_updates',cp.get('material_calls',0))
        return r['execution']['status']=='COMPLETED' and cp.get('completed') and cp.get('cell_count')==1 and coupling>0 and cp.get('checkpoints_loaded')==0 and abs(cp.get('strain_end',-1)-c['protocol']['end_strain'])<1e-7
    finished=[r for r in formal if completed(r)];workflow='PASS' if len(finished)==4 else ('PARTIAL' if formal else 'PENDING')
    optional=lambda name:read_json(base/name) if (base/name).exists() else None
    material=optional('fluid_measured.json');mapping=optional('membrane_measured.json');frozen=optional('formal_frozen.json')
    model=mapping.get('model_comparability','PARTIAL') if mapping else 'PARTIAL'
    screenchecks={};feedback={};prepared=[];sensitivity={};repeats={}
    for engine,emptyname in [('hemocell','empty_hemocell'),('mirheo','empty_dpd')]:
        empty=lookup.get(emptyname,{});fluid=empty.get('fluid')
        if fluid and empty.get('completion',{}).get('completed'):
            gates={name:fluid.get(field) for name,field in [('empty_velocity_relative_l2','profile_relative_l2'),('effective_shear_relative_error','effective_shear_relative_error'),('near_wall_density_relative_error','near_wall_density_relative_error'),('mass_relative_drift','mass_relative_drift')]}
            if 'temperature_relative_error' in fluid:gates['temperature_relative_error']=fluid['temperature_relative_error']
            screenchecks[emptyname]=dict(status='PASS' if all(v<=c['screen'][k] for k,v in gates.items()) and fluid.get('wall_crossings_max',0)==0 else 'FAILED',metrics=gates,details=fluid)
        else:screenchecks[emptyname]=dict(status='NOT_MEASURED')
        primary=lookup.get('main_'+engine+'_1',{});strict=lookup.get('strict_'+engine,{})
        sensitivity[engine]=strict_comparison(primary,strict,c,frozen['strict_comparison_window'] if frozen else [.5,2.])
        screenchecks['strict_'+engine]=sensitivity[engine]
        flow_candidates=[lookup[n] for n in ('main_'+engine+'_1','main_'+engine+'_2') if n in lookup and (lookup[n].get('observed_last_saved_strain') or 0)>=2.]
        flow_source=max(flow_candidates,key=lambda r:r['observed_last_saved_strain']) if flow_candidates else primary
        window=[2.,min(4.,flow_source.get('observed_last_saved_strain') or 0)]
        if window[1]>window[0] and empty.get('directory'):
            cell_local=flow_source.get('local_flow');empty_local=empty.get('local_flow')
            if not cell_local or not np.allclose(cell_local['window'],window,atol=1e-8):cell_local=local_flow_summary(flow_source['directory'],read_csv,window)
            if not empty_local or not np.allclose(empty_local['window'],window,atol=1e-8):empty_local=local_flow_summary(empty['directory'],read_csv,window)
            feedback[engine]=dict(feedback_comparison(cell_local,empty_local,optional('geometry.json')['a']),task=flow_source.get('task'),reference_task=emptyname,window=window)
        else:feedback[engine]=dict(status='NOT_MEASURED')
        rows=[r for r in finished if r['task'].startswith('main_'+engine)]
        values=[r['costs']['total_s'] for r in rows]
        repeats[engine]=dict(count=len(values),independent_cold_processes=True,task_names=[r['task'] for r in rows],
                            end_to_end_s=values,mean_s=float(np.mean(values)) if values else None,
                            minimum_s=min(values) if values else None,maximum_s=max(values) if values else None,
                            sample_std_s=float(np.std(values,ddof=1)) if len(values)>1 else None,
                            seconds_per_strain=[v/c['protocol']['end_strain'] for v in values],
                            coupled_only_s=[r['completion']['coupled_s'] for r in rows],
                            dynamics=[dynamics_window(r,[2.,4.]) for r in rows],confidence_interval='NOT_ESTIMATED: only two cold repetitions planned; deterministic LBM repeats do not sample model uncertainty')
    if material:
        screenchecks['viscosity_half_dt']=dict(status='PASS' if material['strict_difference']<=c['screen']['viscosity_strict_relative_difference'] else 'FAILED',relative_difference=material['strict_difference'],tolerance=c['screen']['viscosity_strict_relative_difference'])
    else:screenchecks['viscosity_half_dt']=dict(status='NOT_MEASURED')
    if mapping:
        screenchecks['native_membrane_response']=dict(status='PASS' if mapping['max_response_relative_difference']<=c['screen']['membrane_response_relative_difference'] else 'FAILED',relative_difference=mapping['max_response_relative_difference'],tolerance=c['screen']['membrane_response_relative_difference'])
    else:screenchecks['native_membrane_response']=dict(status='NOT_MEASURED')
    for r in formal:
        screenchecks[r['task']+'_geometry']=dict(status=r.get('checks',{}).get('geometry_screen','NOT_MEASURED'),details=r.get('checks',{}))
        if r.get('leakage_probes',{}).get('mismatches',0):screenchecks[r['task']+'_fluid_membership']=dict(status='FAILED',details=r['leakage_probes'])
    for repeat in (1,2):
        left=lookup.get('main_hemocell_'+str(repeat),{});right=lookup.get('main_mirheo_'+str(repeat),{})
        initial=lambda r:next((f for f in r.get('frames',[]) if f['phase']=='shear' and abs(f['strain'])<1e-9),None)
        l,r=initial(left),initial(right)
        if l and r:
            x,y=np.array(l['vertices']),np.array(r['vertices']);xc,yc=x.mean(axis=0),y.mean(axis=0);a=optional('geometry.json')['a']
            rms=float(np.sqrt(np.mean(np.sum(((x-xc)-(y-yc))**2,axis=1)))/a)
            prepared.append(dict(repeat=repeat,relative_shape_rms=rms,center_difference=(yc-xc).tolist(),HemoCell=l['metrics'],Mirheo=r['metrics'],
                                 status='PASS' if rms<=c['screen']['initial_relaxed_shape_relative_rms'] else 'FAILED',definition='Corresponding vertex RMS after removing center translation only, normalized by reference a; no rotation or time fitting.'))
    screenchecks['prepared_shapes']=dict(status=('PASS' if all(p['status']=='PASS' for p in prepared) else 'FAILED') if len(prepared)==2 else 'NOT_MEASURED',pairs=prepared)
    runtime_full=lookup.get('diagnostic_mirheo_full',{});runtime_half=lookup.get('diagnostic_mirheo_halfdt',{})
    runtime_complete=bool(runtime_full and completed(runtime_full))
    runtime_diagnostics=dict(plan=optional('runtime_candidate_frozen.json'),completed_full_cold_runs=int(runtime_complete),required_cold_repetitions=2,
                             full_task=runtime_full.get('task'),costs=runtime_full.get('costs'),geometry=runtime_full.get('checks'),leakage=runtime_full.get('leakage_probes'),
                             last_saved_strain=runtime_full.get('observed_last_saved_strain'),execution_status=runtime_full.get('execution',{}).get('status'),
                             sensitivity=strict_comparison(runtime_full,runtime_half,c,[.5,1.]),
                             counts_as_original_repeats=False,status='ONE_FULL_COLD_WORKFLOW' if runtime_complete else 'NOT_COMPLETED_TO_GAMMA4')
    if runtime_full.get('local_flow') and lookup.get('empty_dpd',{}).get('local_flow'):
        runtime_diagnostics['feedback']=dict(feedback_comparison(runtime_full['local_flow'],lookup['empty_dpd']['local_flow'],optional('geometry.json')['a']),
                                            task='diagnostic_mirheo_full',reference_task='empty_dpd',limitation='Same physical inputs, different CUDA queue setting; difference contains thermal scatter and is not a noise-free causal force measurement.')
    screen='FAILED' if any(x['status']=='FAILED' for x in screenchecks.values()) else ('PASS' if all(x['status']=='PASS' for x in screenchecks.values()) and workflow=='PASS' else 'NOT_VERIFIED')
    qualified=speedup(repeats['hemocell']['mean_s'] or 0,repeats['mirheo']['mean_s'] or 0,workflow=workflow,model=model,screen=screen,repeats=min(x['count'] for x in repeats.values()))
    reasons=[]
    if workflow!='PASS':reasons.append('HemoCell 两次主运行到 Γ=4；Mirheo 原配置最后保存到 Γ=0.6、2.9，半步到 Γ=0.7；受限队列诊断最后到 Γ=3.4。均因原生粗碰撞候选溢出退出，未取得 Mirheo 完整 Γ=4 的总时间。' if optional('native_failure_final.json') else '双方两次完整冷启动尚未全部完成。')
    else:reasons.append('双方均完成两次真实双向耦合冷启动到 Γ=4，可以报告各自实际工作流成本。')
    if model!='PASS':reasons.append('独立膜响应、弯曲参考和节点惯性尚未对齐，不能作同等物理问题的性能排名。')
    if screen=='FAILED':reasons.append('事前筛选存在未通过项，未调整阈值或拟合最终轨迹。')
    return dict(schema_version=2,workflow=workflow,model_comparability=model,benchmark_screen=screen,qualified_speedup=qualified,
                human_review='PENDING',research_applicability='NOT_VALIDATED',config=c,geometry=optional('geometry.json'),units=__import__('py_scripts.single_rbc_benchmark.physics',fromlist=['units']).units(),
                common_endpoint_strain=c['protocol']['end_strain'],formal_cold_completed=len(finished),runs=all_runs,
                solver_charged_s=sum(r['execution']['elapsed_monotonic_s'] for r in all_runs if '/solver/' in r['directory']),
                preparation_cpu_s=sum(r['execution']['elapsed_monotonic_s'] for r in all_runs if '/preparation/' in r['directory']),
                one_time_solver_checks_s=sum(r['execution']['elapsed_monotonic_s'] for r in all_runs if '/solver/' in r['directory'] and r not in formal),
                material=material,material_audit=optional('material_diagnostics.json'),membrane_mapping=mapping,dimensionless=optional('dimensionless_audit.json') or optional('dimensionless_measured.json'),formal_frozen=frozen,
                screen_checks=screenchecks,feedback=feedback,repeat_statistics=repeats,prepared_shapes=prepared,sensitivity=sensitivity,runtime_diagnostics=runtime_diagnostics,
                environment=optional('environment_authorized_start_corrected.json') or optional('environment_authorized_start.json'),authorization=optional('authorization.json'),protection=optional('protection_current.json'),cpu_tests=optional('cpu_tests_final.json'),native_failure=optional('native_failure_final.json'),
                historical_install_time='UNKNOWN_FOR_THIS_MATCHED_CASE; software already installed, no reinstall cost invented',
                strict_impermeability='NOT_VERIFIED',space_convergence='NOT_VERIFIED',period_estimate='NOT_ESTABLISHED; no period inferred from incomplete or nonperiodic trajectories',
                reason=''.join(reasons))
