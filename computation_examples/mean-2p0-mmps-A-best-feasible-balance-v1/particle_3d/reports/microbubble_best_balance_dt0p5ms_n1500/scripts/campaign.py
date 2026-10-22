"""Isolated new-flow cohort; reuse protected source and exact finite-MB integrator.

No original module globals are changed. Explicit private bindings change only
the authorized flow identity, environment and nominal time step.
"""
from pathlib import Path
from functools import partial
from types import MethodType
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse, cProfile, gzip, json, multiprocessing, os, sys, time, traceback
import numpy as np

HERE = Path(__file__).resolve().parents[1]
FLOW_SHA = 'fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2'
DT = .0005
N = 1500
JOB = None


def dump(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.pending')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    os.replace(temp, path)


def setup():
    config = json.loads((HERE/'config.json').read_text())
    root = Path(config['source_root'])
    sys.path.insert(0, str(root/'particle_3d/src'))
    from particle_3d.formal_cohort_p9a5 import digest
    from particle_3d.population_inlet_p9a4 import load_new_environment
    from particle_3d.field import FrozenFEMField
    from particle_3d.inlet_flux import frozen_boundary_flux
    from particle_3d.injection_admission import FiniteSizeAdmission
    import pyvista as pv
    snapshot = root/'particle_3d/reports/network_derived_flow_mb_validation_v1/data/code_snapshot_manifest.json'
    protected = json.loads(snapshot.read_text())
    for relative, expected in protected.items():
        assert digest(root/relative) == expected, relative
    flowpath = Path(config['flow_path'])
    assert digest(flowpath) == FLOW_SHA
    env = load_new_environment(root)
    flow = pv.read(flowpath)
    # validate_pair proves canonical geometry, connectivity and numbering agree.
    env.field = FrozenFEMField.from_grids(env.mesh, flow)
    env.audit, samplers = frozen_boundary_flux(env.mesh, flow, env.boundaries)
    env.sampler = samplers['INLET']
    env.checker = FiniteSizeAdmission(None, None, wall=env.wall, field=env.field)
    env.new_flow_path, env.new_flow_sha256 = flowpath, FLOW_SHA
    assert env.mu == .00345312
    return root, env, protected


def prepare():
    root, env, protected = setup()
    from particle_3d.continuous_infusion import ContinuousInfusionSource, PopulationLedger
    from particle_3d.formal_cohort_p9a5 import content_sha, digest, write_new
    from particle_3d.particle6_stepper import bind_query_dependency
    from particle_3d.population_inlet_p9a4 import generate
    original = json.loads((root/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json').read_text())
    contract = dict(schema='BEST_BALANCE_MB1500_DT0P5MS_V1', flow_sha256=FLOW_SHA,
        dt_s=DT, count=N, horizons_s=[3.,6.,12.], master_seed=2026092805,
        source='FIRST_1500_ACCEPTED_IN_SOURCE_ORDER_WITHOUT_OUTCOME_SELECTION',
        source_distribution=original['source_distribution'],
        diameter_min_m=original['D_min_m'], diameter_max_m=original['D_max_m'],
        concentration_m3=original['C_modeled_D_le_4um_m3'],
        concentration_semantics=original['concentration_semantics'],
        concentration_measured=False, Q_in_m3_s=env.sampler.Q_m3_s,
        lambda_source_s_inv=original['C_modeled_D_le_4um_m3']*env.sampler.Q_m3_s,
        old_source_contract_sha256=digest(root/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json'),
        finite_size_admission='ORIGINAL_SINGLE_SHOT_NO_RETRY',
        dynamics='UNCHANGED_P9A1_P65_SINGLE_BUBBLE_ONE_WAY_NO_MB_MB_NO_RBC',
        coverage='NO_OUTCOME_SELECTION_OR_FORCE; REPORT_ACTUAL_OUTLETS',
        authorization='User requested NEW visualized flow, 0.5 ms, about 1500 trajectories',
        GPU_role='FP64 mesh/gradient and trajectory diagnostics; CPU finite-size dynamics',
        CPU_workers=8, protected_source_count=len(protected))
    write_new(HERE/'data/contract.json', contract)
    # A separately named real-data adapter accepts only the verified new hash.
    # It does NOT invoke the synthetic-data bypass or alter the original class.
    class NewFlowSource(ContinuousInfusionSource):
        __init__ = bind_query_dependency(ContinuousInfusionSource.__init__,
                                        {'NEW_FLOW_SHA256': FLOW_SHA})
    source = NewFlowSource(sampler=env.sampler, distribution=env.distribution,
        checker=env.checker, concentration_m3=contract['concentration_m3'],
        flow_sha256=FLOW_SHA, master_seed=contract['master_seed'],
        source_contract_sha256=content_sha(contract), scientific_identity={
            'protected_snapshot_sha256':digest(root/'particle_3d/reports/network_derived_flow_mb_validation_v1/data/code_snapshot_manifest.json'),
            'adapter_sha256':digest(__file__)})
    ledger = PopulationLedger(source.identity)
    while ledger.accepted < N:
        ledger = generate(source, max(512, int((N-ledger.accepted)*5.3)),
                          workers=6, ledger=ledger)
        print('SOURCE', len(ledger.rows), ledger.accepted, flush=True)
    events = ledger.birth_events(N)
    # Only the predeclared accepted prefix is the scientific cohort. Retain
    # extra proposals from the last parallel chunk as unused source records.
    write_new(HERE/'data/cohort.json', dict(count=N, events=events,
        identity=source.identity, selection=contract['source']))
    with gzip.open(HERE/'data/proposal_ledger.jsonl.gz','wt') as stream:
        for row in ledger.rows: stream.write(json.dumps(row, allow_nan=False)+'\n')
    write_new(HERE/'data/preparation.json', dict(flow_sha256=FLOW_SHA,
        flux=env.audit, mu_pa_s=env.mu, source_proposals=len(ledger.rows),
        accepted_proposals=ledger.accepted, selected=N,
        selected_source_event_last=events[-1]['source_event_id'],
        selected_birth_time_last_s=events[-1]['birth_time_s'],
        protected_sources=protected, source_identity=source.identity))
    np.savez_compressed(HERE/'data/gpu_mesh_input.npz', points=env.field.points_m,
        tetra=env.field.tetra, velocity=env.field.velocity_nodes_m_s,
        pressure=env.field.pressure_nodes_pa, inverse=env.field.geometry.inverse,
        gradients=env.field.gradients_s_inv, wall_triangles=env.wall.triangles)
    print('PREPARED', N, flush=True)


def bind_job(env, identity):
    from particle_3d import formal_dynamics_p9a5 as dynamics
    from particle_3d.formal_cohort_p9a5 import HorizonSchedule
    from particle_3d.particle6_stepper import bind_query_dependency
    audit = bind_query_dependency(dynamics.audit_result, {'ENV':env})
    return bind_query_dependency(dynamics.job, dict(ENV=env, IDENTITY=identity,
        DT=DT, audit_result=audit, HorizonSchedule=partial(HorizonSchedule, dt=DT)))


def run_job(spec):
    return JOB(spec)


def run(stage):
    global JOB
    root, env, protected = setup()
    from particle_3d.formal_cohort_p9a5 import digest, completion_matches
    cohort = json.loads((HERE/'data/cohort.json').read_text())
    identity = dict(flow_sha256=FLOW_SHA, dt_s=DT,
        cohort_sha256=digest(HERE/'data/cohort.json'),
        contract_sha256=digest(HERE/'data/contract.json'),
        campaign_sha256=digest(__file__), protected_sources=protected,
        finite_size_integrator='ORIGINAL_P9A1_P65')
    JOB = bind_job(env, identity)
    events = cohort['events']
    if stage == 'pilot':
        profiler = cProfile.Profile(); profiler.enable()
        row = JOB((events[0], str(HERE/'pilot/reference/mb_000001')))
        profiler.disable(); profiler.dump_stats(str(HERE/'logs/pilot.prof'))
        dump(HERE/'data/pilot.json', row); print('PILOT',json.dumps(row),flush=True)
        return
    start=time.time(); rows=[]; errors=[]
    (HERE/'tracks').mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=8, mp_context=multiprocessing.get_context('fork')) as pool:
        futures={pool.submit(run_job,(e,str(HERE/'tracks'/f"mb_{e['particle_id']:06d}"))):e for e in events}
        for f in as_completed(futures):
            event=futures[f]
            try: rows.append(f.result())
            except Exception: errors.append(dict(particle_id=event['particle_id'],error=traceback.format_exc()))
            counts={k:sum(r['outlet']==k for r in rows) for k in ('O1','O2','O3')}
            progress=dict(finished=len(rows),failed_execution=len(errors),total=N,
                outlets=counts,elapsed_s=time.time()-start,
                statuses={s:sum(r['status']==s for r in rows) for s in sorted({r['status'] for r in rows})})
            dump(HERE/'data/progress.json',progress)
            print(json.dumps(progress),flush=True)
    rows.sort(key=lambda r:r['particle_id'])
    dump(HERE/'data/metrics.json',rows);dump(HERE/'data/errors.json',errors)
    verified = all(completion_matches(HERE/'tracks'/f"mb_{e['particle_id']:06d}",identity,e) for e in events)
    safety=['penetration_count','handoff_violation_count','inlet_escape_count','nan_inf_count','unclassified_corruption_count']
    totals={k:sum(r[k] for r in rows) for k in safety}
    summary=dict(progress,all_completion_hashes_verified=verified,safety_totals=totals,
        integration_finished=len(rows)==N and not errors,
        numerical_gate=verified and not any(totals.values()) and not any(r['status']=='SOLVER_FAILURE' for r in rows),
        dt_s=DT,flow_sha256=FLOW_SHA,identity=identity)
    dump(HERE/'data/production_complete.json',summary)
    print('PRODUCTION_COMPLETE',json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','pilot','production']);a=p.parse_args()
    (prepare if a.stage=='prepare' else lambda:run(a.stage))()
