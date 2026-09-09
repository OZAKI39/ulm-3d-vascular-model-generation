"""One predeclared diagnostic, charged to the unchanged shared authorization."""
import copy
import json
from pathlib import Path
import shlex
import shutil
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT,read_json,write_json,sha256_file,atomic_state,environment_identity
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.runner import exclusive_lock,shared_budget_state,run_attempt
from py_scripts.fluid_physics.analysis import read_csv
from py_scripts.fluid_comparison.experiments import load_config as load_comparison
from .calculations import profile_arrays,recover_profile_moments
from .observations import phase_audit,velocity_statistics


def completion_state(execution,completion,required_steps):
    return 'COMPLETED_FIXED_PHYSICAL_WINDOW' if (execution['status']=='COMPLETED' and not execution['timeout'] and
        completion['steps']==required_steps and completion['status']=='COMPLETED_PLANNED_STEPS') else 'PARTIAL_OR_FAILED_NOT_A_COMPLETE_COMPARISON'


def verify_probe_cache(directory,c):
    d=Path(directory)
    for name,h in read_json(d/'output_sha256.json').items():
        p=(d/name).resolve()
        if not p.is_relative_to(d.resolve()) or sha256_file(p)!=h:raise ValueError('PROBE_RAW_HASH_MISMATCH')
    preflight=read_json(d/'preflight.json');package=Path(preflight['diagnostic_package']);verify_package(package)
    if sha256_file(package/'package_sha256.json')!=preflight['diagnostic_manifest_sha256']:raise ValueError('PROBE_PREFLIGHT_IDENTITY_CHANGED')
    provenance=read_json(package/'provenance.json')
    if provenance['identity']['config_sha256']!=c['_config_sha256']:raise ValueError('PROBE_CONFIG_CHANGED_NO_AUTOMATIC_REUSE')
    current=environment_identity()
    for key in ('library_sha256','mirheo_commit','precision'):
        if current[key]!=provenance['environment'][key]:raise ValueError('PROBE_NATIVE_ENVIRONMENT_CHANGED')


def execute(c,package):
    verify_package(package);s=read_json(package/'diagnosis_summary.json');plan=read_json(package/'proposed_experiments.json')['selected_probe']
    cc=load_comparison(c['comparison_config']);campaign=Path(c['runs_root'])/c['campaign_id'];task_id=c['probe']['id']
    worker=Path(__file__).with_name('gpu_worker.py');helper=Path(__file__).with_name('observations.py')
    existing=campaign/task_id
    if existing.exists():
        # A view/replay must never retry even a failed/partial attempt.
        verify_probe_cache(existing,c)
        if sha256_file(existing/'gpu_worker.py')!=sha256_file(worker) or sha256_file(existing/'observations.py')!=sha256_file(helper):
            raise ValueError('PROBE_WORKER_CHANGED_NO_AUTOMATIC_RETRY')
        print('PROBE_CACHE_REUSED_NO_GPU '+str(existing),flush=True);return existing
    validation=read_json(Path(c['setup_directory'])/'CPU_validation.json')
    if validation['status']!='PASS':raise ValueError('CPU_VALIDATION_REQUIRED')
    for path,h in validation['code_sha256'].items():
        if sha256_file(path)!=h:raise ValueError('CPU_VALIDATION_STALE '+path)
    if plan['status']!='ELIGIBLE_WITHIN_EXISTING_AUTHORIZATION':
        print('REQUIRES_ADDITIONAL_AUTHORIZATION: no probe launched',flush=True);return None
    if s['original_reproduction']!='PASS':raise ValueError('REPRODUCTION_REQUIRED')
    env=environment_identity();prior=read_json(package/'provenance.json')['environment']
    for key in ('library_sha256','mirheo_commit','precision'):
        if env[key]!=prior[key]:raise ValueError('NATIVE_ENVIRONMENT_CHANGED')
    old=next(r for r in read_json(package/'corrected_tasks.json') if r['task_id']=='sdpd_equilibrium')
    spec=copy.deepcopy(old['parameters']);dt=plan['dt_star']
    spec.update(dt_star=dt,steps=plan['steps'],desired_steps=plan['steps'],snapshot_every=round(.0002/dt),native_stats_every=round(.0002/dt),
        script_sha256=sha256_file(worker),observation_helper_sha256=sha256_file(helper),prepared_package=str(package),
        source_category='DIAGNOSTIC_TRANSIENT_ONLY',diagnostic_plan=plan)
    spec['task'].update(id=task_id,test_kind='diagnostic_thermal_trend',dt_factor=.5,allocation_s=plan['allocation_s'])
    runner_config=copy.deepcopy(cc);runner_config['campaign_id']=c['campaign_id'];runner_config['require_shared_budget']=True
    pool_path=Path(cc['shared_budget_pool']);pool=read_json(pool_path)
    with exclusive_lock(pool['gpu_lock']):
        state=shared_budget_state(Path(cc['runs_root'])/cc['campaign_id'],cc)
        if state['remaining_s']<plan['allocation_s'] or any(m['running_reservation_s'] for m in state['members']):
            raise RuntimeError('SHARED_BALANCE_CHANGED_NO_LAUNCH')
        pool=state['pool'];campaign.mkdir(parents=True,exist_ok=True)
        if not any(Path(m['directory']).resolve()==campaign.resolve() for m in pool['members']):
            write_json(campaign/'shared_pool_before_registration.json',pool)
            pool['members'].append({'directory':str(campaign),'campaign_id':c['campaign_id'],'ledger_required':False})
            atomic_state(pool_path,pool)
    key={'config_sha256':c['_config_sha256'],'worker':sha256_file(worker),'helper':sha256_file(helper),
         'library':env['library_sha256'],'spec':spec}
    def command(d):
        write_json(d/'actual_parameters.json',spec);write_json(d/'preflight.json',{'plan':plan,'CPU_validation':validation,
          'diagnostic_package':str(package),'diagnostic_manifest_sha256':sha256_file(package/'package_sha256.json'),
          'shared_budget_before_reservation':{k:v for k,v in state.items() if k!='pool'}})
        shutil.copyfile(worker,d/'gpu_worker.py');shutil.copyfile(helper,d/'observations.py');(d/'control').mkdir()
        shutil.copyfile(PROJECT_ROOT/'vendor/Mirheo/LICENSE',d/'MIRHEO_LICENSE.txt')
        launch='#!/usr/bin/env bash\nset -euo pipefail\nsource '+shlex.quote(str(PROJECT_ROOT/'scripts/activate_mirheo.sh'))+'\nexec /usr/bin/mpirun.openmpi --bind-to none -np 2 '+shlex.quote(env['python'])+' -B -u '+shlex.quote(str(d/'gpu_worker.py'))+' --spec '+shlex.quote(str(d/'actual_parameters.json'))+'\n'
        (d/'launch.sh').write_text(launch);return ['/bin/bash',str(d/'launch.sh')]
    print(f"PROBE {task_id} fixed_steps={spec['steps']} dt={dt} allocation={plan['allocation_s']} s",flush=True)
    d,record,hit=run_attempt(campaign,runner_config,task_id,key,command,plan['allocation_s'])
    print('PROBE_RESULT '+json.dumps({'directory':str(d),'status':record['status'],'charged_s':record.get('elapsed_monotonic_s'),'cache':hit}),flush=True)
    return d


def analyze_probe(directory,c,old_rows,u):
    d=Path(directory)
    verify_probe_cache(d,c)
    execution=read_json(d/'execution.json');spec=read_json(d/'actual_parameters.json')
    completion=read_json(d/'worker_completion.json') if (d/'worker_completion.json').exists() else {'steps':0,'status':'MISSING_COMPLETION'}
    status=completion_state(execution,completion,spec['steps'])
    if not (d/'moments.csv').exists() or completion['steps']==0:
        return {'summary':{'status':status,'charged_s':execution['elapsed_monotonic_s'],'interpretation':'No complete physical comparison; no automatic retry.'},'data':None}
    base=next(r for r in old_rows if r['task_id']=='sdpd_equilibrium');bm=read_csv(Path(base['directory'])/'moments.csv')
    m,p,n,v,v2=profile_arrays(d,spec);mom=recover_profile_moments(n,v,v2,spec['m_star'])
    windows=[]
    for lo,hi in c['probe']['predeclared_comparison_windows_star']:
        entry={'interval_star':[lo,hi]}
        for label,record in [('baseline_dt',bm),('probe_half_dt',m)]:
            mask=(record['time_star']>lo)&(record['time_star']<=hi+1e-12);x=record['kBT_COM_star'][mask]
            entry[label]={'samples':len(x),'mean':float(x.mean()) if len(x) else None,'min':float(x.min()) if len(x) else None,'max':float(x.max()) if len(x) else None,
              'complete_interval':bool(record['time_star'][-1]>=hi-1e-12),'CI95':None,'CI_reason':'Single drifting trajectory; intervals are not independent repeats.'}
        windows.append(entry)
    mask=bm['time_star']<=c['probe']['duration_star']+1e-12;bx=bm['kBT_COM_star'][mask];px=m['kBT_COM_star']
    observations=read_json(d/'phase_observations.json') if (d/'phase_observations.json').exists() else []
    cpu=[]
    for i,r in enumerate(observations):
        prefix=d/f'snapshot_{i:02d}'
        arrays=[np.load(str(prefix)+'_'+name+'.npy',allow_pickle=False) for name in ('positions','pre_velocities','pre_forces','post_positions','post_velocities')]
        cpu.append({'step':r['step'],**phase_audit(*arrays,spec['dt_star'],spec['m_star'],spec['domain_star']),
                    'velocity':velocity_statistics(arrays[-1],spec['m_star'])})
    summary={'status':status,'charged_s':execution['elapsed_monotonic_s'],'steps':completion['steps'],
      'physical_duration_s':completion['steps']*spec['dt_star']*u.t0,'baseline_peak':float(bx.max()),'half_dt_peak':float(px.max()),
      'peak_relative_difference':float(px.max()/bx.max()-1),'window_means':windows,
      'moment_COM_max_abs_error':float(np.max(abs(mom['COM']-px))),
      'interpretation':'Both matched unforced runs retain a large early heating peak and decay. Halving dt does not remove that transient; initialization/structure relaxation is supported. Stable late temperature bias is not tested.' if status=='COMPLETED_FIXED_PHYSICAL_WINDOW' and px.max()>2 and bx.max()>2 and windows[-1]['probe_half_dt']['mean']<windows[0]['probe_half_dt']['mean'] else 'Incomplete or different trend: inspect recorded windows; no steady-state conclusion.',
      'qualification_capable':False,'same_noise_path':False,'direct_old_initial_state_hash_available':False,
      'remaining_own_processes':execution['remaining_own_group_processes']}
    density=read_csv(d/'density_samples.csv')
    from .native import snapshot_audit
    structure=[snapshot_audit(d/'initial_positions_global.npy',spec,initial_state=True)]
    frames=sorted(d.glob('snapshot_*_positions.npy'))
    # Exclude separately saved post-kick coordinates: their density is not at
    # the same phase and must never be substituted into this comparison.
    frames=[p for p in frames if '_post_' not in p.name]
    for path in [frames[0],frames[-1]] if frames else []:structure.append(snapshot_audit(path,spec))
    return {'summary':summary,'initial_state':read_json(d/'initial_state.json'),'observations':observations,'CPU_phase_recheck':cpu,
      'structure_audit':structure,
      'channel_layouts':read_json(d/'channel_layouts.json') if (d/'channel_layouts.json').exists() else None,
      'series':{'baseline_time_star':bm['time_star'][mask].tolist(),'baseline_temperature':bx.tolist(),
                'half_dt_time_star':m['time_star'].tolist(),'half_dt_temperature':px.tolist(),
                'density_time_star':density['density_time_star'].tolist(),'kernel_mean':density['kernel_number_mean'].tolist()},
      'raw_directory':str(d),'raw_manifest_sha256':sha256_file(d/'output_sha256.json')}
