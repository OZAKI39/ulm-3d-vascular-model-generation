from collections import Counter
import numpy as np
import pytest
from particle_3d.particle81_simulation import environment,lock_upstream
from particle_3d.particle81_replay import OUTLETS,frame_schedule
from particle_3d.particle8_replay import read,digest,canonical_hash


def test_upstream_and_complete_mesh_preserved(root):
    assert lock_upstream()>1700
    env=environment();geo=np.load(root/'data/full_frozen_geometry.npz')
    for role,mesh in env.boundaries.items():
        np.testing.assert_array_equal(geo[role+'_points_m'],mesh.points)
        np.testing.assert_array_equal(geo[role+'_faces'],mesh.faces.reshape(-1,4)[:,1:])
    assert len(geo['WALL_faces'])==45221


def test_required_dataset_scale_and_terminal_accounting(scene,root):
    cat=scene.catalog;b=read(root/'data/birth_ledger.json')
    assert cat['completed']>=300
    assert set(scene.entries)=={b['particle_id'] for b in b['events']}
    assert sum(cat['outlet_counts'].values())==cat['completed']
    assert len(scene.entries)==sum(cat['end_reasons'].values())
    assert all(e['end_reason'] for e in scene.entries.values())


def test_every_trajectory_raw_samples_and_geometry(scene,root):
    env=environment();budget=env.wall.roundoff_m
    for pid,e in scene.entries.items():
        a=scene.arrays[pid];m=read(root/e['metadata_path'])
        assert digest(root/e['metadata_path'])==e['metadata_sha256']
        assert np.isfinite(a).all() and np.all(np.diff(a[:,0])>0)
        assert len(a)==e['sample_count']
        if len(a):
            np.testing.assert_array_equal(a[0,10:14],e['initial_q'])
            np.testing.assert_array_equal(a[0,1:4],e['inlet_position_m'])
            np.testing.assert_allclose(np.diff(a[:,1:4],axis=0),np.diff(a[:,0])[:,None]*a[1:,4:7],atol=budget,rtol=0)
            assert a[:,15].min()>=-budget*1.01
            assert np.allclose(np.linalg.norm(a[:,10:14],axis=1),1,rtol=0,atol=1e-14)
            assert m['radius_m']==m['birth_metadata']['radius_m']
            assert abs(np.linalg.norm(np.diff(a[:,1:4],axis=0),axis=1).sum()-e['path_length_m'])<1e-18
        else:assert not e['completed'] and e['end_reason']=='INLET_ADMISSION_GUARD_EXHAUSTED_UNRESOLVED'


def test_every_completed_path_reaches_actual_outlet(scene):
    env=environment()
    for pid in scene.completed:
        e=scene.entries[pid];a=scene.arrays[pid]
        assert len(a)>=2 and e['exit_outlet'] in OUTLETS and e['birth_time_s']<e['exit_time_s']
        hit=env.classifier.first_event(a[-2,1:4],a[-1,1:4])
        assert hit is not None and hit.role==e['exit_outlet']
        assert e['path_length_m']>50e-6 # Actual full-domain branch distance, not a 1 ms entrance segment.
        assert np.linalg.norm(a[-1,1:4]-a[0,1:4])>35e-6


def test_birth_proposals_on_official_inlet_and_no_radius_resampling(root,scene):
    tri=np.array(read(root/'data/inlet_mesh.json')['triangles_m'])
    for pid,e in scene.entries.items():
        m=read(root/e['metadata_path'])
        for c in m['admission_candidates']:
            t=tri[c['triangle_id']];p=np.array(c['position_m']);uv=np.linalg.lstsq((t[1:]-t[0]).T,p-t[0],rcond=None)[0]
            assert np.linalg.norm(p-(t[0]+uv@(t[1:]-t[0])))<1e-17
            assert min(*uv,1-uv.sum())>-1e-10
        if e['sample_count']:
            np.testing.assert_allclose(np.array(e['inlet_barycentric'])@tri[e['inlet_face']],e['inlet_position_m'],atol=1e-18,rtol=0)


def test_event_ledger_order_and_accounting(root,scene):
    events=read(root/'data/particle8_1_events.json')['events']
    assert events==sorted(events,key=lambda e:(e['time_s'],e['particle_id'],e['order']))
    for pid,e in scene.entries.items():
        selected=[x for x in events if x['particle_id']==pid]
        assert selected[0]['kind']=='SCHEDULED_BIRTH'
        assert selected[-1]['kind']==('DELETED_AFTER_EXIT' if e['completed'] else 'EXPLICIT_SAFETY_STOP' if e['sample_count'] else 'UNRESOLVED_ADMISSION_TERMINAL_RECORD')
    end=scene.snapshot(scene.catalog['replay_end_time_s']);c=end['counts']
    assert c['active']==0 and c['scheduled']==c['completed']+c['stopped']+c['unresolved']


@pytest.mark.parametrize('kind',['01_full_vessel_replay','02_accumulation','03_single_journeys','04_outlet_ensemble'])
def test_display_timing_identity_and_no_future_history(scene,kind):
    schedule=frame_schedule(scene,kind);times=[r['time_s'] for r in schedule]
    assert np.all(np.diff(times)>=0)
    for spec in schedule:
        state=scene.snapshot(spec['time_s'],spec['only_ids'])
        assert not state['pending_drawn_ids']
        for pid in state['completed_ids']:assert scene.entries[pid]['exit_time_s']<=spec['time_s']
        for p in state['active']:
            e=scene.entries[p['particle_id']]
            assert e['birth_time_s']<=spec['time_s']<e['last_physical_time_s']
            position,q=scene.position(p['particle_id'],p['elapsed_time_s'])
            np.testing.assert_array_equal(position,p['position_m']);np.testing.assert_array_equal(q,p['q'])


def test_interpolation_endpoints_and_no_extrapolation(scene):
    pid=scene.completed[0];a=scene.arrays[pid]
    for k in [0,len(a)//2,len(a)-1]:
        x,q=scene.position(pid,a[k,0]);np.testing.assert_array_equal(x,a[k,1:4]);np.testing.assert_array_equal(q,a[k,10:14])
    for t in [-1e-9,a[-1,0]+1e-9]:
        with pytest.raises(ValueError):scene.position(pid,t)


def test_independent_population_and_physics_audit(root,scene):
    audit=read(root/'data/physics_audit.json')
    assert audit['all_pass'] and audit['all_population_exact']
    assert set(audit['frozen_original_input_sha256_checks'])=={'mesh','flow','WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03'}
    assert all(audit['frozen_original_input_sha256_checks'].values())
    assert audit['source_scene_sha256']==scene.sha256
    assert len(audit['original_population_checks'])==scene.catalog['scheduled']
    for row in audit['original_population_checks']:
        assert all(row[k] for k in ['original_size_exact','original_isotropic_q_exact','deterministic_birth_exact','first_flux_proposal_exact'])
    assert len(audit['completed_final_outlet_checks'])==scene.catalog['completed']
    assert all(p['official_final_segment_match'] for p in audit['completed_final_outlet_checks'])
    assert all(p['inside_original_fem'] for p in audit['field_probes'])
    assert digest(root/audit['derived_sample_export'])==audit['derived_export_sha256']
    assert audit['derived_sample_count']==scene.catalog['physical_samples']
    assert not audit['failures_are_physiological_capture_probabilities']


def test_full_sample_csv_has_original_values_speed_and_clock(root,scene):
    import csv,gzip
    counts=Counter()
    with gzip.open(root/'data/trajectory_samples.csv.gz','rt',newline='') as f:
        for row in csv.DictReader(f):
            pid=int(row['particle_id']);index=counts[pid];counts[pid]+=1
            a=scene.arrays[pid][index];entry=scene.entries[pid]
            for k,col in enumerate(scene.catalog['original_sample_columns']):assert float(row[col])==a[k]
            assert float(row['physical_time_s'])==entry['birth_time_s']+a[0]
            assert np.isclose(float(row['speed_m_s']),np.linalg.norm(a[4:7]),rtol=4*np.finfo(float).eps,atol=0)
            assert int(row['active_until_terminal'])==(index<len(scene.arrays[pid])-1)
    assert sum(counts.values())==scene.catalog['physical_samples']
    for pid,a in scene.arrays.items():assert counts[pid]==len(a)


def test_limited_timestep_sensitivity_keeps_original_births(root):
    audit=read(root/'data/timestep_sensitivity.json')
    assert not audit['production_convergence_established']
    assert audit['physical_cohort_count_excludes_these_reintegrations']
    assert {(r['particle_id'],r['dt_factor']) for r in audit['rows']}=={(pid,f) for pid in [1,4,7,13] for f in [2,4]}
    for r in audit['rows']:
        assert r['same_birth_metadata'] and r['same_inlet'] and r['all_finite']
        m=read(root/r['refined_metadata_path'])
        path=(root/r['refined_metadata_path']).parent.parent/m['samples_path']
        assert digest(path)==r['refined_samples_sha256']
        assert m['end_reason']==r['refined_end']


@pytest.mark.parametrize('target',['birth','geometry','lock','samples','metadata'])
def test_replay_refuses_stale_or_corrupted_saved_inputs(root,scene,tmp_path,target):
    import shutil,json
    from particle_3d.particle81_replay import Scene
    cat=dict(scene.catalog);cat['entries']=[scene.entries[scene.completed[0]]];entry=cat['entries'][0]
    paths={'birth':'data/birth_ledger.json','geometry':'data/full_frozen_geometry.npz',
           'lock':'data/upstream_lock.json','samples':entry['samples_path'],'metadata':entry['metadata_path']}
    for relative in paths.values():
        dest=tmp_path/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/relative,dest)
    (tmp_path/'data/particle8_1_trajectory_catalog.json').write_text(json.dumps(cat))
    # A byte change that leaves most formats parseable must still invalidate the
    # source hash, rather than silently produce a different displayed scene.
    with (tmp_path/paths[target]).open('ab') as f:f.write(b' ')
    with pytest.raises(ValueError,match='changed'):Scene(tmp_path)
