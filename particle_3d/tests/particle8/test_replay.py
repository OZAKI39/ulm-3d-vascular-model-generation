from copy import deepcopy
import math
import numpy as np
import pytest
from particle_3d.particle8_replay import (REPO,P7,read,digest,validate_scene,snapshot,interpolate,slerp,build_bundle)


def test_upstream_original_science_and_evidence_unchanged(output):
    p=read(output/'data/upstream_provenance.json')
    for mapping in ['upstream_reports_sha256','consumed_data_sha256','particle7_source_sha256']:
        for path,h in p[mapping].items(): assert digest(REPO/path)==h,path
    assert p['source_files_verified']==285


@pytest.mark.parametrize('key',['synthetic','restarted','real_mixed','real_single_mb'])
def test_schema_and_frozen_assumptions(scenes,key):
    scene=scenes[key]; assert validate_scene(scene)
    assert scene['assumptions']['H_D']==.45 and scene['assumptions']['C_MB']==8.5e12
    assert scene['assumptions']['cumulative_flux_deterministic']
    assert scene['assumptions']['isotropic_orientation_v0']
    assert 'REAL_RBC_CONTINUOUS_ADMISSION_AND_PASSAGE_NOT_ESTABLISHED' in scene['carried_forward_limitations']


@pytest.mark.parametrize('key,filename',[('synthetic','10_continuous_final.json'),('restarted','10_restarted_final.json'),('real_mixed','13_real_smoke.json')])
def test_all_original_geometry_times_and_provenance_are_preserved(scenes,key,filename):
    original=read(P7/'data'/filename);scene=scenes[key]
    assert len(original['events'])==len(scene['records'])
    for event,r in zip(original['events'],scene['records']):
        assert r['geometry']['original_sample']==event
        assert r['scheduled_time_s']==event['scheduled_time_s']
        assert r['geometry']['orientation_wxyz']==event['q']
        assert r['particle_id']==event['particle_id']
        if event['species']=='RBC': assert r['geometry']['aspect_ratio']==event['geometry']['c_m']/event['geometry']['a_m']


def test_actual_lammps_full_restart_state_exact(scenes,output):
    a=read(P7/'data/10_continuous_final.json');b=read(P7/'data/10_restarted_final.json')
    for k in ['events','births','exits','pending','active','scheduler']: assert a[k]==b[k]
    assert read(output/'data/restart_exact_fields.json')['actual_destroy_and_binary_restore']
    for t in [0,.1249,.125,.1251,.25]: assert snapshot(scenes['synthetic'],t)==snapshot(scenes['restarted'],t)
    assert snapshot(scenes['restarted'],.125)['restart_segment']=='post_restart'


def test_every_exit_deletes_once_no_loss(scenes):
    scene=scenes['synthetic'];ledger=scene['event_ledger']
    ids=[r['particle_id'] for r in scene['records']];assert len(ids)==len(set(ids))==2358
    for r in scene['records']:
        if r['exit_time_s'] is not None:
            t=r['exit_time_s'];assert r['delete_time_s']==t
            before=snapshot(scene,float(np.nextafter(t,-math.inf)))
            after=snapshot(scene,t)
            assert r['particle_id'] in {p['particle_id'] for p in before['active']}
            assert r['particle_id'] not in {p['particle_id'] for p in after['active']}
    assert len([e for e in ledger if e['status']=='deleted'])==2354
    s=snapshot(scene,.25);assert s['counts']['RBC']['active']==4
    assert s['counts']['MB']['deleted']==2 and s['counts']['RBC']['deleted']==2352


def test_synthetic_interpolation_matches_original_common_plug(scenes):
    for r in scenes['synthetic']['records']:
        a,b=r['trajectory'];mid=(a['time_s']+b['time_s'])/2;s=interpolate(r,mid)
        expected=np.array(a['position_m'])+np.array([0,0,2.5e-5])*(mid-a['time_s'])
        np.testing.assert_allclose(s['position_m'],expected,atol=2e-21,rtol=0)
        assert s['q']==a['q']
        assert interpolate(r,a['time_s'])['position_m']==a['position_m']
        assert interpolate(r,b['time_s'])['position_m']==b['position_m']


def test_real_smoke_static_and_pending_is_not_in_lumen(scenes):
    scene=scenes['real_mixed'];s=snapshot(scene,scene['time_range_s'][1])
    assert s['counts']['RBC']==dict(scheduled=1110,admitted=1,pending=1109,active=1,exited=0,deleted=0)
    assert s['counts']['MB']==dict(scheduled=1,admitted=0,pending=1,active=0,exited=0,deleted=0)
    for p in s['pending']: assert p['position_m'] is None
    admitted=scene['records'][0];a,b=admitted['trajectory']
    assert a['position_m']==b['position_m']
    assert interpolate(admitted,(a['time_s']+b['time_s'])/2)['position_m']==a['position_m']
    assert admitted['admitted_shape']['mode']=='CAPILLARY_DEFORMED'
    assert len([r for r in scene['records'] if r['trajectory']])==1


def test_real_single_mb_uses_original_accepted_p65_states(scenes):
    d=read(P7/'data/13_isolated_mb_transport.json');scene=scenes['real_single_mb'];r=scene['records'][0]
    assert r['geometry']['radius_m']==d['geometry']['radius_m']
    for old,new in zip(d['states'],r['trajectory']):
        assert old['time_s']==new['time_s'] and old['particle']['position']==new['position_m']
    delta=np.linalg.norm(np.array(r['trajectory'][-1]['position_m'])-r['trajectory'][0]['position_m'])
    assert delta==d['displacement_m']
    assert not snapshot(scene,scene['time_range_s'][0])['active']
    assert snapshot(scene,r['admit_time_s'])['active'][0]['particle_id']==1


def test_pending_identity_and_wait_duration(output):
    d=read(output/'data/pending_identity.json');a=d['original'];b=d['admitted']
    assert d['identity_exact']
    for k in ['particle_id','geometry','q','provenance','scheduled_time_s']:
        assert a[k]==b[k] and all(r[k]==a[k] for r in d['blocked'])
    assert d['replay_record']['pending_duration_s']==b['admitted_time_s']-a['scheduled_time_s']


@pytest.mark.parametrize('bad',['duplicate','class','time','position','quaternion','pending_position'])
def test_bad_replay_rejected(scenes,bad):
    s=deepcopy(scenes['real_single_mb']);r=s['records'][0]
    if bad=='duplicate': s['records'].append(deepcopy(r))
    elif bad=='class': r['source_classification']='synthetic_control'
    elif bad=='time': r['scheduled_time_s']=float('nan')
    elif bad=='position': r['trajectory'][0]['position_m'][0]=float('nan')
    elif bad=='quaternion': r['trajectory'][0]['q']=[0,0,0,0]
    else:r['admit_time_s']=None
    with pytest.raises(ValueError):validate_scene(s)


def test_no_extrapolation_or_pending_motion(scenes):
    s=scenes['real_single_mb'];r=s['records'][0]
    for t in [s['time_range_s'][0]-1,s['time_range_s'][1]+1,float('nan')]:
        with pytest.raises(ValueError):snapshot(s,t)
    with pytest.raises(ValueError):interpolate(r,r['admit_time_s']-1e-5)
    with pytest.raises(ValueError):interpolate(scenes['real_mixed']['records'][-1],42.986)


def test_slerp_antipodal_and_shortest_arc():
    q=np.array([1.,0,0,0]);np.testing.assert_array_equal(slerp(q,-q,.5),q)
    r=slerp(q,np.array([0.,1,0,0]),.5);np.testing.assert_allclose(r,[2**-.5,2**-.5,0,0],atol=1e-15)


def test_reexport_from_cache_is_reproducible(tmp_path,output):
    build_bundle(tmp_path)
    for name in ['synthetic_scene.json','restarted_scene.json','real_single_mb_scene.json','real_mixed_scene.json','pending_identity.json']:
        assert (tmp_path/'data'/name).read_bytes()==(output/'data'/name).read_bytes()
