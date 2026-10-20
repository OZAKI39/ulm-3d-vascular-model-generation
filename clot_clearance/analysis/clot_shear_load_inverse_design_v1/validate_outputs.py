"""Replay original-ID damage and validate delivery, without new physics runs."""
from analyze import *
def main(case_names=('uniform','counter_rotating')):
    ref=load_npz(O/'provenance/reference.npz');assert read(O/'KEY_CHECKS.json')['status']=='PASS'
    baseline=read(O/'BASELINE_MANIFEST.json')
    for path,item in baseline['files'].items():assert sha(path)==item['sha256'],path
    n=int(ref['clot_count']);N=len(ref['X']);b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];result={}
    for name in case_names:
        check_start=time.perf_counter()
        case=O/'cases'/name;summary=read(case/'SUMMARY.json');request=read(case/'REQUEST.json');rows=histories(case)
        assert summary['status']=='COMPLETED' and len(rows)==summary['accepted_steps'] and abs(rows[-1]['time_s']-END)<1e-14
        # Runtime inputs stay frozen. Independently developed postprocessors are
        # not imported by the formal driver and may be completed after launch.
        runtime={'common.py','native.py','run_cases.py','current_surface.py','shear_field_library.py','multifragment_shear_objective.py'}
        for filename in runtime:assert sha(A/filename)==request['source_sha256'][filename],filename
        postprocessing_changes=[filename for filename,want in request['source_sha256'].items() if sha(A/filename)!=want]
        assert not runtime.intersection(postprocessing_changes)
        for filename,want in request['input_sha256'].items():assert sha(O/filename)==want,filename
        assert request['library_sha256']==sha(LIB) and summary['original_onsets_unchanged'] and summary['native_topology_matches_current_native']
        initial=load_npz(case/'keys/initial.npz');assert np.array_equal(initial['d'],cut.astype(float))
        error=max(abs(q['W_J']-(q['U_J']-initial['scalars'][0]+q['K_J']-initial['K']+q['D_J']+q['Uc_J']-initial['Uc']-q['Jc_J'])-q['R_J']) for q in rows)
        assert error<1e-20 and max(q['relative_R'] for q in rows)<=.01
        assert max(abs(q['W_J']-q['W_drive_J']-q['W_drag_J']) for q in rows)<1e-19
        assert all(q['drag_power_W']<=1e-20 and q['W_drag_J']<=1e-20 for q in rows)
        assert all(abs(q['drive_power_W'])<1e-20 for q in rows if q['time_s']>=OFF)
        assert all(q['free_nodes_with_retained_surface']+q['free_nodes_without_resolved_surface']==q['detached_nodes'] for q in rows)
        assert max(q['max_dd'] for q in rows)<=.1+1e-10 and all(q['numerical_damping_work_J']==0 for q in rows)
        keys={};tool=CurrentSurface(ref,name);force_error=0.;states=0
        for p in (case/'keys').glob('*.npz'):
            assert 'failed' not in p.name
            s=load_npz(p);step=int(s['step']);keys.setdefault(step,[]).append(p)
            assert np.array_equal(s['atom_ids'],np.arange(1,N+1)) and np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any()
            assert np.all(s['d'][cut]==1) and not s['d'][bt==3].any()
            oracle=tool.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))['force'];force_error=max(force_error,float(np.max(abs(oracle-s['external_force']))))
            assert np.max(abs(s['force_stage_velocity']-s['velocity']))<1e-10
            assert force_error<1e-15 and not s['external_force'][n:].any() and float(s['failed_tension_max_N'])<1e-20
            assert np.isfinite(s['x']).all() and np.isfinite(s['kappa']).all() and np.all(s['kappa']>=0)
            states+=1
        registry={int(k):v for k,v in read(case/'COMPONENT_REGISTRY.json').items()};count=0;topology_updates=0
        events=gzip.open(case/'damage_events.jsonl.gz','rt');pending=next(events,None);damage=cut.astype(float);cached=graph.topology(ref,damage);lasttime=-1;previous_free=None
        with gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','rt') as f:
            for line in f:
                z=json.loads(line);step=z['step'];assert z['time_s']>lasttime;lasttime=z['time_s'];changed=False
                if pending:
                    e=json.loads(pending)
                    assert e['step']>=step
                    if e['step']==step:
                        ids=np.asarray(e['original_bond_ids'])-1;values=np.asarray(e['d']);assert len(set(ids))==len(ids) and not cut[ids].any() and np.all(values>=damage[ids]) and np.all(values<=1)
                        changed=bool(np.any((values==1)&(damage[ids]<1)));damage[ids]=values;pending=next(events,None)
                if changed:cached=graph.topology(ref,damage);topology_updates+=1
                lab,sizes,att=cached;actual={tuple((np.flatnonzero(lab==i)+1).tolist()) for i in np.flatnonzero(~att)}
                saved={tuple(registry[q['component_id']]['atom_ids']) for q in z['components']};assert actual==saved
                assert z['attached_nodes']+z['detached_nodes']==n and sum(map(len,saved))==z['detached_nodes']
                if step:assert rows[step-1]['detached_nodes']==z['detached_nodes'] and abs(rows[step-1]['time_s']-z['time_s'])<1e-15
                key_state=None
                for p in keys.get(step,[]):
                    state=load_npz(p);assert np.array_equal(state['d'],damage),p
                    if key_state is None:key_state=state
                if changed or previous_free!=saved:
                    for members in saved:
                        boundary,inside=graph.boundary(ref,np.asarray(members)-1);assert not np.any(damage[boundary]<1)
                for q in z['components']:
                    assert not q['anchor_reachable'] and q['remaining_live_external_connections']==0
                    if key_state is not None:
                        s=key_state;ids=np.asarray(registry[q['component_id']]['atom_ids'])-1
                        sh=shape(s['x'][ids]);assert np.allclose(sh['actual_span_m'],q['actual_span_m'],rtol=0,atol=1e-18)
                        _,inside=graph.boundary(ref,ids);inside=inside[(bt[inside]==1)&~cut[inside]]
                        assert q['internal_initial_bonds']==len(inside) and q['internal_alive_bonds']==int(np.sum(damage[inside]<1))
                previous_free=saved;count+=1
        events.close();assert pending is None and count==len(rows)+1
        for e in summary['secondary_events']:
            parent=set(registry[e['parent_component_id']]['atom_ids']);children=[set(registry[q]['atom_ids']) for q in e['children_component_ids']]
            assert set.union(*children)==parent and sum(map(len,children))==len(parent)
        for event in ['first_softening','first_failure']:
            e=summary['events'].get(event)
            if e:
                s=load_npz(case/e['evidence']);ids=np.asarray(e['original_bond_ids'])-1
                assert np.all(s['kappa'][ids]>S0[bt[ids]]) and not cut[ids].any()
                assert np.all(s['d'][ids]==1) if event=='first_failure' else np.all(s['d'][ids]>0)
        evidence=[case/filename for filename in ['history.csv','damage_events.jsonl.gz','CONNECTED_COMPONENT_HISTORY.jsonl.gz','COMPONENT_REGISTRY.json','REQUEST.json','SUMMARY.json']]+sorted((case/'keys').glob('*.npz'))
        result[name]=dict(status='PASS',accepted_records=len(rows),component_records=count,topology_recomputations=topology_updates,accepted_FP64_states=states,all_accepted_memberships_match_original_graph=True,energy_reconstruction_error_J=error,current_surface_force_error_N=force_error,runtime_source_and_protocol_frozen=True,nonruntime_postprocessing_files_completed_after_launch=postprocessing_changes,validator_sha256=sha(A/'validate_outputs.py'),raw_evidence_sha256={str(p.relative_to(case)):sha(p) for p in evidence},validation_wall_s=time.perf_counter()-check_start)
        dump(case/'DELIVERY_VALIDATION.json',result[name])
    combined={}
    for name in ['uniform','counter_rotating']:
        case=O/'cases'/name;p=case/'DELIVERY_VALIDATION.json'
        if not p.exists():continue
        proof=read(p);assert proof['validator_sha256']==sha(A/'validate_outputs.py')
        for filename,want in proof['raw_evidence_sha256'].items():assert sha(case/filename)==want
        combined[name]=proof
    status='PASS' if len(combined)==2 else 'PARTIAL_ONE_CASE'
    dump(O/'DELIVERY_VALIDATION.json',dict(status=status,cases=combined,scope='Original-ID damage replay, full accepted component histories, shape, native shear coupling, accounting, frozen source/protocol. Each case validated once; unchanged proofs combined by SHA. No new physics case.'))
    print('DELIVERY_VALIDATION',status,{k:{q:v[q] for q in ['accepted_records','component_records','accepted_FP64_states','validation_wall_s']} for k,v in combined.items()})
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--case',choices=['uniform','counter_rotating']);args=p.parse_args();main([args.case] if args.case else ['uniform','counter_rotating'])
