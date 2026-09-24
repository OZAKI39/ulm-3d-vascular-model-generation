#!/usr/bin/env python3
"""Gated, isolated P9-A.2 inlet sampling and unchanged P9-A.1 integration."""
import argparse,json,hashlib,time,os,socket,platform
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from collections import Counter
import numpy as np
from particle_3d.injection_method_c import MethodCSource,TruncatedSonoVue,METHOD,LEGACY_METHOD,canonical_bytes,sha
from particle_3d.particle81_simulation import environment,dump,DT,HORIZON
from particle_3d.particle9a_provenance import require_current_flow,REPO
from particle_3d.particle8_replay import canonical_hash
from particle_3d.particle7_cases import SONOVUE

R=REPO/'particle_3d/reports/particle9a2_inlet_sampling';D=R/'data'
O=REPO/'particle_3d/outputs/particle9a2_2mmps'
SOURCE=None;POINT=None
ROLES=['OUTLET_01','OUTLET_02','OUTLET_03','NO_EXIT']

def read(p):return json.loads(Path(p).read_text())

def setup(seed):
    global SOURCE
    e=environment();flow=require_current_flow()
    geometry={k:v['sha256'] for k,v in e.provenance['boundary_manifest']['boundaries'].items() if k in ('INLET','WALL')}
    SOURCE=MethodCSource(sampler=e.sampler,wall=e.wall,field=e.field,distribution=TruncatedSonoVue(SONOVUE),
        seed=seed,flow_sha256=flow['flow_sha256'],geometry_sha256=geometry)
    return SOURCE

def sample_shard(ids):
    return [SOURCE.event(i) for i in ids]

def point_shard(spec):
    events,folder=spec;rows=[]
    for e in events:
        trace=POINT.trace(e['birth_center_m']);path=trace.pop('path')
        assert np.array_equal(path[0,1:],e['birth_center_m'])
        p=Path(folder)/f"point_{e['particle_id']:06d}.npz";p.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(p,path=path)
        rows.append(dict(particle_id=e['particle_id'],point_outlet=trace['outlet'] or 'NO_EXIT',
            point_end_reason=trace['end_reason'],diameter_m=e['diameter_m'],position_m=e['birth_center_m'],
            birth_metadata_sha256=canonical_hash(e),path_sha256=sha(p),sample_count=len(path),
            initial_position_exact=True,terminal_classifier_verified=bool(trace['outlet'])))
    return rows

def source_identity():
    return {p.name:sha(p) for p in sorted((REPO/'particle_3d/src/particle_3d').glob('*.py'))}

def generate(label,count,seed,workers,folder):
    s=setup(seed);start=time.time();events=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for f in as_completed([pool.submit(sample_shard,list(range(a,min(a+32,count+1)))) for a in range(1,count+1,32)]):
            events.extend(f.result());print(label,'sampled',len(events),'/',count,flush=True)
    events.sort(key=lambda e:e['particle_id'])
    ledger=dict(schema='P9A2_METHOD_C_BIRTH_LEDGER_V1',method=METHOD,events=events,seed=seed,
        flow_sha256=s.flow_sha256,geometry_sha256=s.geometry_sha256,source_distribution_contract_sha256=s.distribution.contract_sha256,
        Q_in_m3_s=s.sampler.Q_m3_s,C_MB_m3=8.5e12,birth_rate_semantics='ACTUAL_ENTERING_C_MB_TIMES_Q',
        entering_distribution='SONOVUE_MAX4UM_CONDITIONED_ON_CURRENT_RIGID_SPHERE_INLET_PASSABILITY')
    path=folder/'admission/birth_ledger.json';path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():raise FileExistsError(path)
    path.write_bytes(canonical_bytes(ledger))
    # Reordered full regeneration checks ALL events, with one process and a fresh tree.
    replay=setup(seed);recreated={i:canonical_bytes(replay.event(i)) for i in range(count,0,-1)}
    assert all(recreated[e['particle_id']]==canonical_bytes(e) for e in events)
    report=dict(stage=label,count=count,workers=workers,seed=seed,wall_seconds=time.time()-start,
        global_rejections=sum(e['diameter_global_rejections'] for e in events),
        source_draws=sum(len(e['source_diameter_draws']) for e in events),
        position_proposals=sum(e['position_proposal_count'] for e in events),
        maximum_proposals=max(e['position_proposal_count'] for e in events),
        invalid_births=0,diameter_changes_during_position_retry=0,
        all_events_byte_identical_reversed_single_worker_replay=True,ledger_sha256=sha(path),
        np_version=np.__version__,python=platform.python_version(),host=socket.gethostname())
    for e in events:
        p,status,_=s.checker.check(e,np.array(e['birth_center_m']),{})
        assert p is not None and status=='ACCEPTED'
        assert all(a['diameter_float64_hex']==np.float64(e['diameter_m']).tobytes().hex() for a in e['position_attempts'])
    dump(D/(label+'_sampling.json'),report)
    return ledger

def run_points(label,events,workers,folder):
    global POINT
    from particle_3d.particle82_point_native import NativePointTracer
    POINT=NativePointTracer(environment());start=time.time();rows=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for f in as_completed([pool.submit(point_shard,(events[i:i+32],str(folder))) for i in range(0,len(events),32)]):
            rows.extend(f.result())
    rows.sort(key=lambda r:r['particle_id'])
    result=dict(stage=label,rows=rows,counts={k:sum(r['point_outlet']==k for r in rows) for k in ROLES},
        point_integrator_source_sha256=POINT.source_sha256,workers=workers,wall_seconds=time.time()-start)
    dump(D/(label+'_point.json'),result);print(label,'point',result['counts'],flush=True)

def integrate_job(spec):
    e,folder,identity=spec
    from particle_3d.particle9a_motion import Particle9AStepper
    from particle_3d.particle6_stepper import bind_query_dependency
    from particle_3d.particle82a_integration import integrate_admitted
    env=environment();created=[];start=time.time()
    def factory(*args,**kwargs):
        s=Particle9AStepper(*args,gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv,**kwargs)
        created.append(s);return s
    integrate=bind_query_dependency(integrate_admitted,{'SavedTrajectoryStepper':factory})
    before=canonical_hash(e);m=integrate(e,output=Path(folder))
    assert before==canonical_hash(e)==m['birth_metadata_sha256']
    if not created:raise ValueError('P9A2 requires fresh integration, no old cache')
    m.update(dynamics='P9A1_UNCHANGED_DYNAMICS_WITH_METHOD_C_INLET',p9a2_identity=identity,
        receipt=dict(host=socket.gethostname(),pid=os.getpid(),started_unix_s=start,ended_unix_s=time.time()),
        entry_transition_role='OPEN_INLET_PLANE_BIRTH_NO_STREAMLINE_REPOSITIONING',
        timestep_role='UNCHANGED_P9A1_FORMAL_CONFIGURATION')
    path=Path(folder)/'trajectories'/f"mb_{e['particle_id']:06d}.json";dump(path,m)
    return dict(particle_id=e['particle_id'],completed=m['completed'],outlet=m['exit_outlet'] or 'NO_EXIT',
        reason=m['end_reason'],failure_detail=m['failure_detail'],wall_seconds=m['wall_seconds'])

def trajectories(label,ledger,workers,folder):
    events=ledger['events'];start=time.time();identity=dict(source_sha256=source_identity(),flow=require_current_flow(),
        method=METHOD,dt_s=DT,horizon_s=HORIZON,source_contract_sha256=ledger['source_distribution_contract_sha256'])
    config=dict(stage=label,workers=workers,identity=identity,started_unix_s=start,host=socket.gethostname(),server_root=str(REPO))
    dump(D/(label+'_config.json'),config);rows=[];environment()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for f in as_completed([pool.submit(integrate_job,(e,str(folder),identity)) for e in events]):
            row=f.result();rows.append(row);print(label,len(rows),'/',len(events),json.dumps(row),flush=True)
            dump(D/(label+'_progress.json'),sorted(rows,key=lambda r:r['particle_id']))
    dump(D/(label+'_completed.json'),dict(config,rows=sorted(rows,key=lambda r:r['particle_id']),wall_seconds=time.time()-start))

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['audit','smoke','formal']);p.add_argument('--workers',type=int,default=6);a=p.parse_args()
    if not str(REPO).startswith('/workspace/particle9a2_inlet_'):raise ValueError('Isolated server root required')
    if a.workers!=6:raise ValueError('Audited server worker configuration is 6')
    D.mkdir(parents=True,exist_ok=True)
    if a.stage=='audit':
        s=setup(2026092492)
        contract=REPO/'particle_3d/contracts/SONOVUE_DIAMETER_MAX4UM_V1.json'
        contract.write_bytes(canonical_bytes(s.distribution.contract))
        dump(D/'inlet_size_capacity.json',s.capacity)
        d=s.distribution
        lo,hi=s.capacity['D_geometry_max_bracket_m']
        dump(D/'source_probability.json',dict(original_above_4um=1-d.mass,
            global_impossible_source_probability_bracket=[float(1-d.cdf(hi)),float(1-d.cdf(lo))],
            source_max_um=4.,geometry_max_um=lo*1e6))
        np.savez_compressed(D/'inlet_geometry.npz',triangles_m=s.sampler.triangles,q_m_s=s.sampler.q)
        ledger=generate('audit2000',2000,2026092492,a.workers,R/'audit2000')
        run_points('audit2000',ledger['events'],a.workers,R/'audit2000/point')
        # Probability definition and full reproducibility are tested locally first.
        dump(D/'inlet_audit_gate.json',dict(audit_complete=True,source_contract_sha256=sha(contract),
            sample_ledger_sha256=sha(R/'audit2000/admission/birth_ledger.json'),tests_required_before_smoke=True))
    elif a.stage=='smoke':
        gate=read(D/'local_gate.json')
        assert gate['inlet_unit_tests_pass'] and read(D/'inlet_audit_gate.json')['audit_complete']
        ledger=generate('smoke30',30,2026092493,a.workers,R/'smoke30')
        trajectories('smoke30',ledger,a.workers,R/'smoke30')
    else:
        gate=read(D/'smoke30_gate.json');assert gate['no_new_systematic_failure'] is True
        assert read(D/'independent_flux_reference.json')['all_position_agreement'] is True
        assert read(D/'inlet_audit_gate.json')['audit_complete'] is True
        assert read(D/'local_gate.json')['inlet_unit_tests_pass'] is True
        ledger=generate('formal500',500,2026092494,a.workers,O)
        run_points('formal500',ledger['events'],a.workers,O/'point')
        trajectories('formal500',ledger,a.workers,O)
if __name__=='__main__':main()
