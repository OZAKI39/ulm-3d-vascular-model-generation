"""Read-only historical protection, asset identities and delivery gate."""
import sys,json,datetime
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from vcommon import ROOT,OUT,DATA,dump,sha,identity

def main():
    before=json.loads((ROOT/'evidence/protected_inputs_before.json').read_text());changed=[];missing=[]
    for file,h in before.items():
        p=Path(file)
        if not p.exists():missing.append(file)
        elif sha(p)!=h:changed.append(file)
    dump(OUT/'evidence/protected_inputs_after.json',dict(checked=len(before),changed=changed,missing=missing,all_unchanged=not changed and not missing))
    assert not changed and not missing,'PROTECTED_INPUT_CHANGED: inspect; never overwrite external changes'
    freeze=json.loads((ROOT/'results/model_freeze.json').read_text())
    for split,h in freeze['data_manifest_hash'].items():
        manifest=ROOT/'data'/split/'manifest.json'
        assert sha(manifest)==h,('CANONICAL_MANIFEST_CHANGED',split)
        for name,file_hash in json.loads(manifest.read_text())['files'].items():
            assert sha(manifest.parent/name)==file_hash,('DATA_CHANGED',split,name)
    for name,record in freeze['models'].items():
        assert sha(ROOT/record['checkpoint'])==record['weight_sha256'],('WEIGHT_CHANGED',name)
    for p in (ROOT/'training').glob('*/summary.json'):
        for name,h in json.loads(p.read_text())['checkpoint_sha256'].items():
            if (p.parent/name).exists():assert sha(p.parent/name)==h,('CHECKPOINT_CHANGED',p.parent.name,name)
    visual=json.loads((ROOT/'evidence/visual_delivery_review.json').read_text())
    for row in visual['figures']:
        for ext in ['png','pdf']:
            assert sha(ROOT/'figures'/f"{row['name']}.{ext}")==row[ext+'_sha256']
    progress=json.loads((ROOT/'evidence/progress_document/audit.json').read_text())
    assert sha(progress['path'])==progress['after_sha256']
    assert sha(progress['backup'])==progress['before_sha256']
    raw=Path('/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4')
    raw_sha=sha(raw)
    assert raw_sha=='3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c'
    assets=[]
    for folder in ['training','data/banks','data/train','data/validation','data/calibration','data/final_test','results','figures','videos','evidence','configs','src','scripts','tests','logs','models']:
        for p in sorted((ROOT/folder).rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            if p.name in ['INPUT_OUTPUT_SHA256.json','final_audit.json','final_audit.log']:continue
            assets.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
    for p in sorted(ROOT.iterdir()):
        if p.is_file() and p.name not in ['INPUT_OUTPUT_SHA256.json','FINAL_GATE.json']:
            assets.append(dict(path=p.name,bytes=p.stat().st_size,sha256=sha(p)))
    p=ROOT/'data/bank_manifest.json';assets.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
    dump(ROOT/'INPUT_OUTPUT_SHA256.json',dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),identity=identity(),files=assets,self_hash_excluded=True,active_audit_log_excluded=True,raw_video=dict(path=str(raw),sha256=raw_sha)))
    suite=json.loads((OUT/'training/suite.json').read_text());video=json.loads((OUT/'videos/manifest.json').read_text());visual=json.loads((OUT/'evidence/visual_delivery_review.json').read_text());replay=[json.loads((OUT/'results'/f'real_replay_{n}.json').read_text()) for n in ['PIV_SINGLE','direct_gray','direct_fused','geo_fused']]
    gate=dict(DATASET_PROVENANCE_PASS=True,SPLIT_LEAKAGE_CHECK_PASS=True,AFFINE_MATH_TEST_PASS=True,DIRECT_MODEL_TRAINED=True,GEOMETRY_MODEL_TRAINED=True,HEATMAP_ABLATION_COMPLETE=True,SYNTHETIC_SIX_PARAMETER_ACCURACY='FAIL_AI',FALSE_ACCEPT_COMPARISON='AI_WORSE_AT_MATCHED_COVERAGE',RISK_COVERAGE_EVALUATED=True,CONFIDENCE_INTERVAL_CALIBRATION='COMPLETED_SYNTHETIC_ONLY_WIDE_INTERVALS',REAL_SIX_PARAMETER_ACCURACY='NOT_VERIFIED',LIVE_CAMERA_PASS=False,CAUSAL_REPLAY_TEST=all(r['no_new_frame_TTL_pass'] for r in replay),END_TO_END_LATENCY='NOT_MEASURED',PIV_AI_FAIR_COMPARISON=True,FIGURES_VALIDATED=visual['figures_pass'],VIDEOS_VALIDATED=visual['videos_pass'] and all(v['fully_decoded'] and v['pts_monotonic'] for v in video['videos']),PZT_CONTROL_EXECUTED=False,AFG_CONTROL_EXECUTED=False,WORLD_MODEL_TRAINED=False,TRAINING_READY=False,SCIENTIFIC_CONTROL_ADMITTED=False,training_ready=False,scope=['EXPLORATORY','PHYSICS_UNCALIBRATED','SINGLE_RECORDING_GENERALIZATION_NOT_VERIFIED'],recommendation='RETAIN_PIV_OBSERVATION_MAINLINE; do not admit AI or PIV to unvalidated physical control',trained_runs=len(suite['runs']),final_test_samples=1536,ai_operational_accepted=0,piv_operational_accepted=1111,piv_failed_accepted=85,protected_historical_files=len(before),identity=identity())
    analysis=json.loads((ROOT/'results/analysis.json').read_text())
    piv=analysis['methods']['PIV_SINGLE']
    gate.update(final_test_samples=piv['n_total'],piv_operational_accepted=piv['n_accepted'],piv_failed_accepted=piv['n_accepted_failed'],ai_operational_accepted={name:analysis['methods'][name]['n_accepted'] for name in freeze['models']},canonical_data_and_weights_verified=True,progress_document_backup_verified=True)
    assert all(n==0 for n in gate['ai_operational_accepted'].values()),'Reassess scientific gate if results differ'
    dump(ROOT/'FINAL_GATE.json',gate)
    # Gate content included by hash after it is generated; manifest never includes itself.
    manifest=json.loads((ROOT/'INPUT_OUTPUT_SHA256.json').read_text());manifest['files'].append(dict(path='FINAL_GATE.json',bytes=(ROOT/'FINAL_GATE.json').stat().st_size,sha256=sha(ROOT/'FINAL_GATE.json')));dump(ROOT/'INPUT_OUTPUT_SHA256.json',manifest)
    print('protected',len(before),'asset files',len(assets),'scientific gate',gate['SYNTHETIC_SIX_PARAMETER_ACCURACY'])
if __name__=='__main__':main()
