"""Immutable comparison preparation; reuse the shared monotonic GPU runner."""
from pathlib import Path
import math
import shlex
import shutil
import yaml
from py_scripts.fluid_physics.common import (PROJECT_ROOT,config_read,output,read_json,write_json,
    sha256_file,fingerprint,environment_identity,resource_snapshot,atomic_state,now)
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.runner import exclusive_lock,run_attempt,used_budget,shared_budget_state
from .models import frozen_targets,candidates,native_contract


def load_config(path):
    p=Path(path);p=p if p.is_absolute() else PROJECT_ROOT/p
    c=config_read(p);case,mapping,u=frozen_targets(c)
    old=yaml.safe_load((Path(c['historical_result'])/'fluid_physics.yaml').read_text())
    c['proposed_acceptance']=old['proposed_acceptance']
    c['legacy_windows']=read_json(Path(c['historical_result'])/'legacy_acceptance_mapping.json')
    for key in ['n_star']:
        if c['space'][key]!=mapping[key]:raise ValueError('PRIMARY_DISCRETIZATION_CHANGED')
    c['space'].update(m_star=mapping['m_star'],kBT_star=mapping['kBT_star'])
    c['resolved_candidates']=candidates(c,mapping)
    c['calibration']={'tasks':c['experiments']}
    return c


def verify_historical(c):
    old=Path(c['historical_campaign']);lp=old/'budget_ledger.json'
    if sha256_file(lp)!=c['historical_ledger_sha256']:raise ValueError('HISTORICAL_BUDGET_IDENTITY_CHANGED')
    ledger=read_json(lp)
    if ledger['campaign_limit_s']!=3600 or any(a['status']=='RUNNING' for a in ledger['attempts']):
        raise ValueError('UNVERIFIED_HISTORICAL_BUDGET')
    for a in ledger['attempts']:
        d=Path(a['directory'])
        if not d.resolve().is_relative_to(old.resolve()):raise ValueError('HISTORICAL_RUN_OUTSIDE_CAMPAIGN')
        for name,h in read_json(d/'output_sha256.json').items():
            p=(d/name).resolve()
            if not p.is_relative_to(d.resolve()) or sha256_file(p)!=h:raise ValueError('HISTORICAL_RAW_HASH_MISMATCH '+str(p))
        r=read_json(d/'execution.json')
        if r['elapsed_monotonic_s']!=a['charged_s']:raise ValueError('HISTORICAL_CHARGE_MISMATCH')
    return ledger


def register_budget(c,old):
    campaign=output(Path(c['runs_root'])/c['campaign_id']);pool=output(c['shared_budget_pool'])
    lock=Path(c['historical_campaign']).parent/'.gpu_exclusive.lock'
    with exclusive_lock(lock):
        if pool.exists():
            state=shared_budget_state(campaign,c)
            if state['pool']['historical_ledger_sha256']!=c['historical_ledger_sha256']:raise ValueError('BUDGET_SOURCE_CONFLICT')
            return state
        obj={'schema_version':1,'created_at':now(),'limit_s':old['campaign_limit_s'],'gpu_lock':str(lock),
             'historical_ledger_sha256':c['historical_ledger_sha256'],'historical_initial_used_s':used_budget(old),
             'members':[{'directory':c['historical_campaign'],'campaign_id':old['campaign_id'],'ledger_required':True},
                        {'directory':str(campaign),'campaign_id':c['campaign_id'],'ledger_required':False}],
             'accounting':'Sum actual monotonic charges or full unresolved reservations across every member. Original ledger remains unmodified.'}
        write_json(pool,obj)
        return shared_budget_state(campaign,c)


def prepare(c):
    case,mapping,u=frozen_targets(c);old=verify_historical(c);contract=native_contract();env=environment_identity()
    prior=read_json(Path(c['historical_result'])/'provenance.json')['environment']
    for key in ['library_sha256','mirheo_commit','precision']:
        if env[key]!=prior[key]:raise ValueError('NATIVE_ENVIRONMENT_CHANGED '+key)
    budget=register_budget(c,old)
    critical=[PROJECT_ROOT/'py_scripts/fluid_comparison'/n for n in ['models.py','experiments.py','gpu_worker.py']]
    critical += [PROJECT_ROOT/'py_scripts/fluid_physics'/n for n in ['runner.py','units.py','common.py']]
    identity={'config_sha256':c['_config_sha256'],'historical_result_manifest':c['historical_manifest_sha256'],
              'historical_ledger':c['historical_ledger_sha256'],'environment':env,'native_contract':contract,
              'candidates':c['resolved_candidates'],'code_sha256':{str(p):sha256_file(p) for p in critical}}
    key=fingerprint(identity);d=output(Path(c['data_root'])/('prepared_'+key[:16]))
    if d.exists():verify_package(d);return d
    d.mkdir(parents=True,exist_ok=False)
    budget_public={k:v for k,v in budget.items() if k!='pool'}
    design=[]
    for a in c['resolved_candidates']:
        if a['method']!='SDPD':continue
        cs=a['sound_speed_star'];rho=mapping['rho_star'];nu=mapping['required_nu_star'];spacing=mapping['particle_spacing_m']/u.L0
        pressure_changes={k:u.to_star(p,'pressure')/(rho*cs**2) for k,p in case['outlet_gauge_pressures_pa'].items()}
        design.append({'candidate_id':a['id'],'category':'DERIVED_DESIGN_NOT_MEASURED','input_mu_star':a['viscosity_mu_star'],'target_nu_star':nu,
          'sound_speed_si':u.to_si(cs,'velocity'),'estimated_relative_density_changes':pressure_changes,
          'input_EOS_background_at_global_density_pa':u.to_si(cs**2*(rho-a['rho_0_star']),'pressure'),
          'pressure_reference':'Unknown thermodynamic absolute pressure. Compare each measured equilibrium mechanical reference plus the unchanged legacy gauge differences. Background force effect measured independently.',
          'step_estimates_star':{'acoustic_0.1_spacing_over_cs':.1*spacing/cs,'viscous_0.02_spacing_squared_over_nu':.02*spacing**2/nu,
                                 'thermal_displacement_0.05_spacing_over_3vth':.05*spacing/(3*math.sqrt(a['kBT_star']/a['m_star']))},
          'selected_dt_star':a['dt_star'],'selected_dt_s':a['dt_star']*u.t0,'requires_actual_step_validation':True})
    for name,value in [('physical_case',case),('unit_mapping',mapping),('native_sdpd_contract',contract),('resolved_config',c),
                       ('candidates',c['resolved_candidates']),('design_estimates',design),('initial_shared_budget',budget_public),
                       ('preflight_resources',resource_snapshot()),('provenance',{'identity':identity,'created_at':now()})]:
        write_json(d/(name+'.json'),value)
    shutil.copyfile(c['_config_path'],d/'fluid_model_comparison.yaml')
    write_json(d/'package_sha256.json',{p.name:sha256_file(p) for p in d.iterdir() if p.is_file()});verify_package(d)
    return d


def execute(c,prepared,*,only=None,retry_failed=False):
    verify_package(prepared);env=read_json(Path(prepared)/'provenance.json')['identity']['environment']
    case,mapping,u=frozen_targets(c);campaign=output(Path(c['runs_root'])/c['campaign_id']);rows=[]
    worker=PROJECT_ROOT/'py_scripts/fluid_comparison/gpu_worker.py'
    by={a['id']:a for a in c['resolved_candidates']};speed=None
    plan_ids={t['id'] for t in c['experiments']}
    if only and not set(only)<=plan_ids:raise ValueError('TASK_NOT_IN_FROZEN_PLAN')
    for task in c['experiments']:
        if only and task['id'] not in only:continue
        # A previously completed probe supplies speed without executing its worker.
        probe=campaign/'sdpd_probe/worker_completion.json'
        if speed is None and probe.is_file():
            z=read_json(probe);speed=(z['compute_wall_s']+z['snapshot_wall_s'])/z['steps']
        if task['kind']!='probe' and speed is None:raise RuntimeError('SDPD_PROBE_REQUIRED')
        a=by[task['candidate']];dt=a['dt_star']*task['dt_factor']
        every=max(1,round(c['sampling']['sample_interval_star']/dt))
        state=shared_budget_state(campaign,c);allocation=min(task['allocation_s'],state['remaining_s'])
        desired=task.get('desired_steps',int(task.get('desired_time_star',0)/dt));steps=desired
        if speed and task['kind']!='probe':
            usable=max(0,allocation-18-c['budget']['graceful_margin_s'])
            steps=min(steps,int(usable/(speed*c['sampling']['speed_safety_factor'])),c['sampling']['max_steps'])
        steps=max(every,steps//every*every)
        spec={'schema_version':2,'task':task,'candidate':a,'dt_star':dt,'m_star':a['m_star'],'kBT_star':a['kBT_star'],'rc_star':a['rc_star'],
              'domain_star':c['sampling']['domain_star'],'snapshot_every':every,'native_stats_every':every,'y_bin_star':c['space']['output_y_bin_star'],
              'steps':steps,'desired_steps':desired,'measured_preflight_s_per_step':speed,'initial_velocity_seed':c['sampling']['initial_velocity_seed'],
              'locked_units':u.as_dict(),'precision':env['precision'],'script_sha256':sha256_file(worker),'libmirheo_sha256':env['library_sha256'],
              'source_category':'INPUT','method':a['method'],'prepared_package':str(prepared)}
        key={'config':c['_config_sha256'],'task':task,'candidate':a,'units':u.as_dict(),'script':sha256_file(worker),
             'library':env['library_sha256'],'precision':env['precision'],'sampling':c['sampling']}
        def command(d):
            write_json(d/'actual_parameters.json',spec);shutil.copyfile(worker,d/'gpu_worker.py');(d/'control').mkdir()
            shutil.copyfile(PROJECT_ROOT/'vendor/Mirheo/LICENSE',d/'MIRHEO_LICENSE.txt')
            launch='#!/usr/bin/env bash\nset -euo pipefail\nsource '+shlex.quote(str(PROJECT_ROOT/'scripts/activate_mirheo.sh'))+'\nexec /usr/bin/mpirun.openmpi --bind-to none -np 2 '+shlex.quote(env['python'])+' -B -u '+shlex.quote(str(d/'gpu_worker.py'))+' --spec '+shlex.quote(str(d/'actual_parameters.json'))+'\n'
            (d/'launch.sh').write_text(launch);return ['/bin/bash',str(d/'launch.sh')]
        print(f"TASK {task['id']} method={a['method']} steps={steps} dt*={dt} allocation<={allocation:.3f}s",flush=True)
        d,r,hit=run_attempt(campaign,c,task['id'],key,command,allocation,retry_failed=retry_failed)
        print(f"RESULT {task['id']} cache={hit} status={r['status']} wall_s={r.get('elapsed_monotonic_s')}",flush=True)
        rows.append({'task_id':task['id'],'directory':str(d) if d else None,'cache_hit':hit,'execution':r})
        if d and r.get('remaining_own_group_processes'):raise RuntimeError('OWN_COMPUTE_PROCESSES_REMAIN')
        if r['status'] in ['BUDGET_EXHAUSTED','RESOURCE_RESERVE_BLOCKED','HOST_MEMORY_RESERVE_BLOCKED']:break
        if task['test_kind']=='probe' and r['status']!='COMPLETED':break
    return rows
