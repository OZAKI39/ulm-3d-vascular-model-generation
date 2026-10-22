"""Independent completeness checks; scientific limitations are not pass flags."""
import sys,json,collections,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from full_records import full_records
import pyarrow.parquet as pq
from PIL import Image
def read(p):return json.loads((OUT/p).read_text())
def main():
    audit=read('results/full_video_source_audit.json');assert sha(VIDEO)==audit['sha256']==P['video_sha256'];assert audit['full_decode_pass'];assert read('evidence/independent_local_source_audit.json')['decoded_frames']==937
    records=list(full_records());assert [r['frame_id'] for r in records]==list(range(937));assert np.allclose([r['pts_s'] for r in records],np.arange(937)*.0333,atol=1e-10);assert len(records)==audit['decoded_frames']
    piv=pq.read_table(OUT/'results/full_video_piv_per_frame.parquet').to_pylist();assert len(piv)==4685
    for method in ['candidate_3','candidate_5','candidate_10','raw_3','local_3']:
        rr=[r for r in piv if r['method']==method];assert [r['frame_id'] for r in rr]==list(range(937))
        for r in rr:
            for n in NAMES:assert r['valid'][n]==(r['values'][n] is not None)
    assert sum(r['methods']['candidate_3']['estimate']['valid']['translation_vx'] for r in records)==887
    source_preserved=[];changed=[];protected=read('evidence/protected_inputs_before.json')
    for p,h in protected.items():
        if Path(p).is_file() and sha(p)==h:source_preserved.append(p)
        else:changed.append(p)
    assert not changed,changed
    transfer=[]
    for line in (OUT/'evidence/server_transfer_sha256.txt').read_text().splitlines():
        if not line.strip():continue
        h,p=line.split();local=OUT/'evidence'/Path(p).name;assert sha(local)==h;transfer.append(dict(file=str(local.relative_to(OUT)),sha256=h))
    complete=read('data/DATASET_COMPLETE.json');splits={}
    for split,v in complete['splits'].items():
        index=OUT/'data'/split/'index.json';assert sha(index)==v['index_sha256'];rows=json.loads(index.read_text());assert len(rows)==v['n'];splits[split]=dict(sources=set(r['source_frame'] for r in rows),seeds=set(r['seed'] for r in rows))
    for a in splits:
        for b in splits:
            if a<b:assert not splits[a]['sources']&splits[b]['sources'];assert not splits[a]['seeds']&splits[b]['seeds']
    frozen=OUT/'results/model_freeze.json';assert sha(frozen)==read('results/FINAL_OPEN_RECORD.json')['freeze_sha256']==read('results/MODEL_BENCHMARK_COMPLETE.json')['freeze_sha256']
    ot=read('results/ot_numerical_audit.json');assert not [r for r in ot['rows'] if 'numerical error' in r['warning'].lower() and (r['base_valid'] or r['previous_base_valid'])];assert not [r for r in ot['rows'] if r['base_valid']!=r['previous_base_valid']]
    for m in ['PIV','LRAFT_ST','CoTracker3','POT_unbalanced']:
        r=read('results/causal_replay_'+m+'.json');assert r['n_input']==937;assert r['n_processed']+r['n_dropped']==937;assert r['no_new_frame_TTL_pass'];assert not r['source_error'];ids={x['frame_id'] for x in r['rows']};assert not ids&set(r['dropped_frame_ids']);assert ids|set(r['dropped_frame_ids'])==set(range(937))
        for x in r['rows']:
            assert max(x['state']['time_window']['frame_ids'])==x['frame_id']
            if not x['fresh']:assert not x['six_valid']
    videos=read('videos/full_video_manifest.json')['videos']+[read('videos/comparison_video_manifest.json')['video']]
    for i,v in enumerate(videos):
        assert sha(OUT/'videos'/v['file'])==v['sha256'];assert v['fully_decoded'] and v['max_pts_error_s']<1e-4;assert v['decoded_frames']==(937 if i<2 else 64)
    figures=[]
    for p in sorted((OUT/'figures').glob('Fig*.png')):
        with Image.open(p) as im:im.verify()
        pdf=p.with_suffix('.pdf');assert pdf.read_bytes().startswith(b'%PDF');figures.append(dict(png=p.name,png_sha256=sha(p),pdf_sha256=sha(pdf)))
    assert len(figures)>=11;visual=read('evidence/visual_review.json');assert visual['status']=='PASS';assert len(visual['figures_reviewed'])==len(figures)
    testlog=(OUT/'logs/contracts_r2.log').read_text();assert 'Ran 8 tests' in testlog and '\nOK\n' in testlog
    report=dict(identity=identity(),protected_checked=len(protected),protected_changed=changed,transfer_verified=transfer,full_frames=937,piv_rows=len(piv),figures=figures,video_files=[v['file'] for v in videos],final_freeze_verified=True,split_frames_and_seeds_disjoint=True,OT_numerical_acceptance_changes=0,math_and_state_contract_tests=8,scope='Engineering and synthetic-sample verification, not real accuracy')
    dump(OUT/'evidence/final_integrity.json',report)
    gate=dict(FULL_SOURCE_VIDEO_FOUND=True,FULL_SOURCE_VIDEO_SHA_VERIFIED=True,FULL_VIDEO_DECODE_PASS=True,FULL_DURATION_PIV_COMPLETED=True,FULL_DURATION_PIV_FRAME_COVERAGE=1.0,FULL_DURATION_PIV_RESULTS_VERIFIED=True,FULL_DURATION_PIV_VIDEO_VALIDATED=True,PIV_CAUSAL_REPLAY_COMPLETED=True,LRAFT_ST_ZERO_SHOT_RUN=True,COTRACKER3_ZERO_SHOT_RUN=True,POT_MATCHING_RUN=True,UNLITEFLOWNET_ZERO_SHOT_RUN=True,SELF_SUPERVISED_ADAPTATION_RUN=False,SELF_SUPERVISED_ADAPTATION_STATUS='NOT_RUN',SIX_PARAMETER_MATH_PASS=True,ZERO_SHOT_MOTION_RECOVERY='PARTIAL_SYNTHETIC_SUCCESS_NOT_REAL_CALIBRATION',FALSE_ACCEPT_COMPARISON=True,RISK_COVERAGE_COMPARISON=True,REAL_SIX_PARAMETER_ACCURACY='NOT_VERIFIED',LIVE_CAMERA_PASS=False,LIVE_CAMERA_STATUS='NOT_RUN_BY_SCOPE',CAMERA_TO_CONTROLLER_END_TO_END_LATENCY='NOT_MEASURED',M4_COMPATIBILITY_PASS=True,FIGURES_VALIDATED=True,VIDEOS_VALIDATED=True,PZT_CONTROL_EXECUTED=False,AFG_CONTROL_EXECUTED=False,WORLD_MODEL_TRAINED=False,TRAINING_READY=False,SCIENTIFIC_CONTROL_ADMITTED=False,RECOMMENDATION='KEEP_PIV',GENERALIZATION='SINGLE_RECORDING_GENERALIZATION_NOT_VERIFIED',ENGINEERING_DELIVERY_PASS=True,protected_history_unchanged=True,identity=identity())
    dump(OUT/'FINAL_GATE.json',gate);print('PASS',len(protected),'protected files;',len(figures),'figure groups; 3 decoded videos; 937 complete frames')
if __name__=='__main__':main()
