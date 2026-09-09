"""CPU preparation, locked units/candidates, adaptive small GPU campaign."""
from datetime import datetime, timezone
from pathlib import Path
import math
import shutil
import shlex
import numpy as np
from .common import PROJECT_ROOT, output, read_json, write_json, sha256_file, fingerprint, environment_identity, resource_snapshot, atomic_state, now
from .legacy_case import freeze_case
from .units import Units, thermal_density_units, consistency, viscosity_only_alternative
from .boundary_plan import build_boundary_plan
from .runner import run_attempt, ledger_read, used_budget
from .analysis import analyze_task


def package_units(case,c):
    s=c['space'];u=thermal_density_units(s['L0_m'],case['density_kg_m3'],s['n_star'],s['m_star'],s['kBT_star'],case['temperature_K'])
    spacing=s['L0_m']/s['n_star']**(1/3);rc=s['rc_star']*u.L0;N=case['geometry_volume_m3']*s['n_star']/u.L0**3
    # These are component estimates, not measured allocations or initialized particles.
    shellN=case['wall_area_m2']*rc/spacing**3
    buffers=sum(p['existing_extension']['original_area_um2']*1e-12*p['existing_extension']['actual_axial_length_um']*1e-6 for p in case['ports'])
    mapping={**u.as_dict(),'category':'DERIVED_FROM_LOCKED_DESIGN_AND_SOURCE_DENSITY_AND_USER_TEMPERATURE',
             'm_star':s['m_star'],'n_star':s['n_star'],'rho_star':s['m_star']*s['n_star'],'kBT_star':s['kBT_star'],
             'thermostat_temperature_K':u.temperature_K(s['kBT_star']),
             'required_nu_star':u.to_star(case['kinematic_viscosity_m2_s'],'kinematic_viscosity'),
             'required_mu_star':u.to_star(case['dynamic_viscosity_pa_s'],'dynamic_viscosity'),
             'viscosity_status':'TO_BE_MEASURED_IN_THE_SAME_UNITS; units cannot force material agreement',
             'density_thermal_check':consistency(u,n_star=s['n_star'],m_star=s['m_star'],kBT_star=s['kBT_star'],rho_si=case['density_kg_m3'],temperature_K=case['temperature_K']),
             'particle_spacing_m':spacing,'interaction_cutoff_m':rc,'old_geometry_reference_m':s['old_geometry_reference_m'],
             'output_bin_m':s['output_y_bin_star']*u.L0,'future_sdf_spacing_m':None,
             'bubble_to_spacing_ratio':[d/spacing for d in case['future_constraints']['bubble_diameter_m']],
             'bubble_to_rc_ratio':[d/rc for d in case['future_constraints']['bubble_diameter_m']],
             'known_source_min_diameter_to_spacing':2*case['geometry_scales']['source_centerline_radius_min_um']*1e-6/spacing,
             'min_measured_port_Dh_to_spacing':case['geometry_scales']['measured_port_min_hydraulic_diameter_m']/spacing,
             'resolution_scope':'Pure fluid calibration cannot guarantee 2-4 um object gaps or adhesion resolution. Whole lumen minimum and RBC size unknown.',
             'resource_estimates':{'status':'ESTIMATE_NOT_INITIALIZED_NOT_GPU_PEAK','lumen_volume_including_existing_extensions_m3':case['geometry_volume_m3'],
                  'liquid_particles_including_existing_extensions':N,'existing_extension_volume_estimate_m3':buffers,
                  'extension_particles_subset_not_additional':buffers*s['n_star']/u.L0**3,'wall_shell_thickness_m':rc,'wall_shell_particles_estimate':shellN,
                  'liquid_and_wall_particle_array_bytes_range':[(N+shellN)*b for b in (96,256)],
                  'cell_list_and_halo_bytes_range':[(N+shellN)*b for b in (64,256)],
                  'field_statistics_bytes_per_bin_range':[40,128],
                  'temporary_buffer_estimate':'1-2 times particle + cell-list storage, plus CUDA/MPI runtime; full bounding-box SDF not estimated without domain design',
                  'scale_law':'Halving L0 at fixed n* roughly multiplies liquid count by 8, wall shell count by 4; stored fields depend on independent bin spacing.'}}
    candidates=[]
    for cand in c['candidates']:
        row={**cand,'m_star':s['m_star'],'n_star':s['n_star'],'rho_star':s['m_star']*s['n_star'],
             'rc_star':s['rc_star'],'kBT_star':s['kBT_star'],'sigma_star':math.sqrt(2*cand['gamma']*s['kBT_star']),
             'native_random_force_amplitude_star':math.sqrt(2*cand['gamma']*s['kBT_star']/cand['dt_star']),
             'native_random_convention':'sigma=sqrt(2 gamma kBT/dt); pairwise unit-variance noise and wr=(1-r/rc)^power; dissipation wr^2; do not add another thermostat/random force',
             'dt_s':cand['dt_star']*u.t0,'cutoff_m':rc,'particle_mass_kg':s['m_star']*u.M0,
             'number_density_m3':s['n_star']/u.L0**3,'a_unit':'reduced force','gamma_unit':'reduced mass/time','power_unit':'1',
             'parameter_status':'DESIGN_CHOICE_NOT_VALIDATED','locked_unit_mapping':u.as_dict()}
        candidates.append(row)
    return u,mapping,{'status':'UNVERIFIED','candidates':candidates,'selection_requires':'Simultaneous density, temperature, measured viscosity, dt/force sensitivity and pressure-range realizability'}


def prepare(c):
    case,acceptance,source_hashes,g=freeze_case(c)
    env=environment_identity();resources=resource_snapshot();u,mapping,candidates=package_units(case,c)
    identity={'config_sha256':c['_config_sha256'],'source_sha256':source_hashes,'environment':env,'unit_mapping':u.as_dict()}
    key=fingerprint(identity);campaign=output(Path(c['runs_root'])/c['campaign_id']);campaign.mkdir(parents=True,exist_ok=True)
    ref=campaign/'prepared_package.json'
    supersedes=None
    if ref.exists():
        old=read_json(ref)
        old_directory=Path(old['directory']);verify_package(old_directory)
        if old['identity']==key:return old_directory
        old_identity=read_json(old_directory/'provenance.json')['identity']
        # Supplementary source audit may add immutable provenance, while all
        # numerical inputs and existing source hashes remain fixed. Preserve the
        # old package and ledger; do not alter any GPU attempt/cache identity.
        added_only=(old_identity['config_sha256']==identity['config_sha256'] and old_identity['environment']==identity['environment'] and old_identity['unit_mapping']==identity['unit_mapping']
                    and all(source_hashes.get(p)==h for p,h in old_identity['source_sha256'].items()))
        if not added_only:raise ValueError('PREPARED_INPUT_CHANGED: cannot silently reuse or reset campaign')
        supersedes=str(old_directory)
    directory=output(Path(c['data_root'])/('prepared_'+key[:16]));directory.mkdir(parents=True,exist_ok=False)
    for name,obj in [('physical_case',case),('legacy_acceptance_mapping',acceptance),('unit_mapping',mapping),('dpd_candidates',candidates),
                     ('boundary_plan',build_boundary_plan(case,c,u)),('preflight_resources',resources),
                     ('provenance',{'identity':identity,'identity_sha256':key,'source_sha256':source_hashes,'environment':env,'created_at':now(),
                                    'code_reuse':'Stage-one verified loader, safe paths and hashes; mathematical legacy acceptance definitions; native periodic/virial API wiring; offline Plotly.',
                                    'protected_snapshot':'test_code/outputs/fluid_physics/setup_20260908T154306_330273Z/protected_before.json'})]:
        write_json(directory/(name+'.json'),obj)
    shutil.copyfile(c['_config_path'],directory/'fluid_physics.yaml')
    write_json(directory/'package_sha256.json',{p.name:sha256_file(p) for p in directory.iterdir() if p.is_file()})
    if supersedes:
        atomic_state(ref,{'identity':key,'directory':str(directory),'supplements':supersedes,'reason':'Added read-only historical cut, gauge origin and accepted flux-denominator evidence. Numerical inputs and GPU budget unchanged.'})
    else:write_json(ref,{'identity':key,'directory':str(directory)})
    verify_package(directory)
    return directory


def verify_package(directory):
    directory=Path(directory)
    for name,h in read_json(directory/'package_sha256.json').items():
        p=(directory/name).resolve()
        if not p.is_relative_to(directory.resolve()) or sha256_file(p)!=h:raise ValueError('PACKAGE_HASH_MISMATCH '+name)
        if p.suffix=='.json':read_json(p)


def run_campaign(c,prepared,*,retry_failed=False,only=None):
    prepared=Path(prepared);case=read_json(prepared/'physical_case.json');um=read_json(prepared/'unit_mapping.json');u=Units(um['L0'],um['M0'],um['t0'])
    env=environment_identity();candidate={a['id']:a for a in c['candidates']};worker=PROJECT_ROOT/'py_scripts/fluid_physics/gpu_worker.py'
    campaign=output(Path(c['runs_root'])/c['campaign_id']);speed=None;results=[]
    # Read measured probe speed from the charged existing attempt when resuming.
    ledger=ledger_read(campaign,c['campaign_id'],c['budget'])
    for a in ledger['attempts']:
        p=Path(a['directory'])/'worker_completion.json'
        if p.is_file():
            z=read_json(p);speed=(z['compute_wall_s']+z['snapshot_wall_s'])/max(z['steps'],1)
            break
    for task in c['calibration']['tasks']:
        if only and task['id'] not in only:continue
        if task['kind']!='probe' and speed is None:
            print('PREFLIGHT_SPEED_UNAVAILABLE: refusing unmeasured long task',flush=True);break
        cand=candidate[task['candidate']];dt=cand['dt_star']*task['dt_factor'];cal=c['calibration'];s=c['space']
        desired=task.get('desired_steps',int(task.get('desired_time_star',0)/dt))
        steps=desired
        if speed:
            usable=max(0,task['wall_allocation_s']-30-c['budget']['graceful_margin_s'])
            steps=min(desired,int(usable/(speed*cal['speed_safety_factor'])),cal['max_steps'])
        steps=max(cal['snapshot_every'],steps//cal['snapshot_every']*cal['snapshot_every'])
        spec={'schema_version':1,'task':task,'candidate':cand,'dt_star':dt,'m_star':s['m_star'],'kBT_star':s['kBT_star'],'rc_star':s['rc_star'],
              'domain_star':cal['domain_star'],'snapshot_every':cal['snapshot_every'],'native_stats_every':cal['native_stats_every'],
              'y_bin_star':s['output_y_bin_star'],'steps':steps,'desired_steps':desired,'measured_preflight_s_per_step':speed,
              'initial_velocity_seed':cal['initial_velocity_seed'],'locked_units':u.as_dict(),'precision':env['precision'],
              'script_sha256':sha256_file(worker),'libmirheo_sha256':env['library_sha256'],
              'source_category':'DESIGN_CHOICE_EXECUTED; results recorded separately as MEASURED'}
        cache_identity={'config':c['_config_sha256'],'task':task,'candidate':cand,'units':u.as_dict(),
                        'script':sha256_file(worker),'library':env['library_sha256'],'precision':env['precision'],
                        'domain':cal['domain_star'],'sampling':cal,'environment_sha256':env['sha256']}
        def command(d):
            write_json(d/'actual_parameters.json',spec);shutil.copyfile(worker,d/'gpu_worker.py');(d/'control').mkdir()
            shutil.copyfile(PROJECT_ROOT/'vendor/Mirheo/LICENSE',d/'MIRHEO_LICENSE.txt')
            script='#!/usr/bin/env bash\nset -euo pipefail\nsource '+shlex.quote(str(PROJECT_ROOT/'scripts/activate_mirheo.sh'))+'\nexec /usr/bin/mpirun.openmpi --bind-to none -np 2 '+shlex.quote(env['python'])+' -B -u '+shlex.quote(str(d/'gpu_worker.py'))+' --spec '+shlex.quote(str(d/'actual_parameters.json'))+'\n'
            with (d/'launch.sh').open('x') as f:f.write(script)
            return ['/bin/bash',str(d/'launch.sh')]
        print(f"TASK {task['id']} steps={steps} dt*={dt} allocation<={task['wall_allocation_s']}s",flush=True)
        d,rec,cached=run_attempt(campaign,c,task['id'],cache_identity,command,task['wall_allocation_s'],retry_failed=retry_failed)
        print(f"RESULT {task['id']} cached={cached} {rec.get('status')} elapsed={rec.get('elapsed_monotonic_s')}s",flush=True)
        if d:
            results.append(str(d));completion=d/'worker_completion.json'
            if task['kind']=='probe' and completion.exists():
                z=read_json(completion);speed=(z['compute_wall_s']+z['snapshot_wall_s'])/max(z['steps'],1)
            if rec['remaining_own_group_processes']:raise RuntimeError('OWN_PROCESS_GROUP_NOT_EMPTY')
        if rec['status'] in ('BUDGET_EXHAUSTED','RESOURCE_RESERVE_BLOCKED','HOST_MEMORY_RESERVE_BLOCKED'):break
    return results
