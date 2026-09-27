"""Post-run sensitivity, resume and immutable-input evidence, without new cohorts."""
from pathlib import Path
from collections import Counter
import argparse,csv,json,socket,time
import numpy as np
from .particle82_provenance import atomic_json,sha256,require_remote
from .particle82a_analysis import descriptive


def read(path):return json.loads(Path(path).read_text())


def complete_audit(root,provenance,admission_provenance,trajectory_provenance):
    root=Path(root);report=root/'review';start=time.time()
    require_remote(provenance,provenance['hostname'])
    trajectory=root/'trajectories';cohorts=read(trajectory/'TRAJECTORY_COHORTS.json')['cohorts']
    sensitivity={}
    for factor in [2,4]:
        rows=[]
        for method in 'ABC':
            for event in cohorts[method][:12]:
                name=f'mb_{event["particle_id"]:06d}'
                base=trajectory/'formal'/method/'trajectories'/name
                fine=trajectory/'sensitivity'/f'dt_div_{factor}'/method/'trajectories'/name
                b=read(base.with_suffix('.json'));f=read(fine.with_suffix('.json'))
                bs=np.load(base.with_suffix('.npz'))['samples'];fs=np.load(fine.with_suffix('.npz'))['samples']
                assert np.array_equal(bs[0,1:4],fs[0,1:4])
                common_end=min(bs[-1,0],fs[-1,0]);times=bs[bs[:,0]<=common_end,0]
                bi=bs[bs[:,0]<=common_end,1:4]
                fi=np.column_stack([np.interp(times,fs[:,0],fs[:,k]) for k in [1,2,3]])
                delta=np.linalg.norm(bi-fi,axis=1)*1e6
                rows.append(dict(method=method,event_id=event['particle_id'],base_outlet=b['exit_outlet'],
                    refined_outlet=f['exit_outlet'],base_end_reason=b['end_reason'],refined_end_reason=f['end_reason'],
                    same_outlet=b['exit_outlet']==f['exit_outlet'],same_end_reason=b['end_reason']==f['end_reason'],
                    residence_difference_s=(f['residence_time_s']-b['residence_time_s']) if b['completed'] and f['completed'] else None,
                    common_time_max_position_difference_um=float(delta.max()),
                    interpretation='POSITION_COMPARED_ONLY_OVER_SHARED_SAVED_TIME_INTERVAL'))
        sensitivity[str(factor)]=dict(rows=rows,compared=len(rows),outlet_changes=sum(not r['same_outlet'] for r in rows),
            end_reason_changes=sum(not r['same_end_reason'] for r in rows),
            common_time_max_position_difference_um=descriptive([r['common_time_max_position_difference_um'] for r in rows]))
    atomic_json(report/'TIMESTEP_SENSITIVITY.json',dict(variants=sensitivity,
        selection='FIRST_12_ADMITTED_IDS_PER_METHOD_PREDECLARED_NO_OUTCOME_SELECTION',
        base_dt_s=.00025,REMOTE_SERVER_COMPUTE=True,production_results_unchanged=True))
    sweep=read(report/'METHOD_C_ENTRY_SWEEP_AUDIT.json');entry={r['event_id']:r for r in sweep['rows']}
    contingency={}
    with (report/'common_inlet_audit_ledger.csv').open() as f:
        for event in csv.DictReader(f):
            basin=event['basin'];current=event['center_on_plane_accepted']=='True'
            found=int(event['event_id']) in entry
            certified=found and entry[int(event['event_id'])]['handoff_clear']
            row=contingency.setdefault(basin,Counter())
            row['scheduled']+=1;row['current_on_plane_accepted']+=int(current)
            row['C_birth_found']+=int(found);row['C_entry_certified']+=int(certified)
            row['current_rejected_but_C_entry_certified']+=int(not current and certified)
            row['current_accepted_but_C_no_birth']+=int(current and not found)
            row['current_accepted_but_C_entry_not_certified']+=int(current and not certified)
    atomic_json(report/'PAIRED_CURRENT_VS_C_ENTRY.json',dict(by_anchor_basin=contingency,
        comparison='SAME_EVENT_SAME_FIRST_SIZE_SAME_ANCHOR; CURRENT_FIRST_TRIAL_VS_C_FULL_SWEPT_PATH',
        separate_from_method_A_retry_acceptance=True,REMOTE_SERVER_COMPUTE=True))
    groups={}
    for flag in [True,False]:
        selected=[e for e in cohorts['C'] if entry[e['particle_id']]['handoff_clear']==flag]
        records=[read(trajectory/'formal/C/trajectories'/f'mb_{e["particle_id"]:06d}.json') for e in selected]
        groups[str(flag)]=dict(count=len(records),completed=sum(r['completed'] for r in records),
            outlet_counts={b:sum(r['exit_outlet']==b for r in records) for b in ['OUTLET_01','OUTLET_02','OUTLET_03']},
            end_reasons=dict(Counter(r['end_reason'] for r in records)))
    atomic_json(report/'METHOD_C_OUTCOMES_BY_ENTRY_CERTIFICATE.json',dict(groups=groups,
        diagnostic_stratification_only=True,no_births_or_outcomes_removed=True,
        geometry_certificate_not_a_validation_of_finite_size_entry_dynamics=True))
    # Invoke resume using the original source provenance. No full-run receipt is rewritten.
    from . import particle82a_pipeline as admission_runner
    admission_runner.initialize(root/'admission',admission_provenance,read(root/'admission/RUN_CONFIG.json'))
    afiles=sorted((root/'admission/events').glob('events_0000001_0000064.*'))
    before={str(p):dict(sha256=sha256(p),mtime_ns=p.stat().st_mtime_ns) for p in afiles}
    admission_runner.event_shard((1,64,'events'))
    assert all(v==dict(sha256=sha256(Path(p)),mtime_ns=Path(p).stat().st_mtime_ns) for p,v in before.items())
    from . import particle82a_trajectories as runner
    runner.initialize(trajectory,trajectory_provenance);trajectory_resume=[]
    for method in 'ABC':
        event=cohorts[method][0];files=sorted((trajectory/'formal'/method/'trajectories').glob(f'mb_{event["particle_id"]:06d}.*'))
        original={str(p):dict(sha256=sha256(p),mtime_ns=p.stat().st_mtime_ns) for p in files}
        runner.worker((method,event,'formal',.00025))
        assert all(v==dict(sha256=sha256(Path(p)),mtime_ns=Path(p).stat().st_mtime_ns) for p,v in original.items())
        trajectory_resume.append(dict(method=method,event_id=event['particle_id'],files=original))
    atomic_json(report/'RESUME_VERIFICATION.json',dict(all_pass=True,admission_shard=before,
        trajectories=trajectory_resume,verified='EXISTING_FILES_REUSED_WITH_IDENTICAL_BYTES_AND_MTIMES',
        REMOTE_SERVER_COMPUTE=True,hostname=socket.gethostname()))
    snapshots={}
    for name,prov in [('admission',admission_provenance),('trajectories',trajectory_provenance),('completion_audit',provenance)]:
        repo=Path(prov['_repo']);checked={}
        for group in ['source_sha256','frozen_input_sha256','external_scientific_input_sha256']:
            hashes=prov.get(group,{})
            mismatches=[k for k,v in hashes.items() if sha256(repo/k)!=v]
            assert not mismatches,(name,group,mismatches)
            checked[group]=dict(count=len(hashes),mismatches=mismatches)
        snapshots[name]=dict(commit=prov['source_git_commit'],checks=checked)
    atomic_json(report/'REMOTE_IMMUTABILITY_VERIFICATION.json',dict(all_pass=True,snapshots=snapshots,
        REMOTE_SERVER_COMPUTE=True,hostname=socket.gethostname(),seconds=time.time()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True)
    for name in ['provenance','admission-provenance','trajectory-provenance']:p.add_argument('--'+name,required=True)
    a=p.parse_args();provs=[]
    for value in [a.provenance,a.admission_provenance,a.trajectory_provenance]:
        prov=read(value);prov['_repo']=str(Path(value).parent/'repo');provs.append(prov)
    complete_audit(a.root,*provs)
