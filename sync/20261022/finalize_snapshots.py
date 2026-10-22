#!/usr/bin/env python3
from export_snapshots import *

def load_export(name):
    e=Export(name);s=e.root/'sync/20261022'
    for key,file in [('files','COPIED_FILES'),('excluded','EXCLUDED_FILES'),('pruned','PRUNED_DIRECTORIES'),('aliases','SYMLINK_POLICY'),('identities','SOURCE_IDENTITIES')]:
        setattr(e,key,json.loads((s/(file+'.json')).read_text()))
    return e

def compact(e):
    kept=[]
    for r in e.files:
        p=e.root/r['path'];s=r['path'];why=None
        if e.root.name=='simulation':
            if '/f10ag/series/' in s and p.suffix=='.npz':why='derived_plot_series; canonical_F10AG_fixed_ID_dataset_analysis_code_and_final_figures_retained'
            if '/motion_solver_r3/review/all_active_inputs/' in s and p.suffix=='.npz':why='repeated_legacy_review_arrays; comparison_metrics_failure_evidence_and_original_runs_retained'
        if why:
            e.excluded.append(dict(source=r['source'],path=r['path'],source_bytes=r['source_bytes'],sha256=r['source_sha256'],reason=why));p.unlink();continue
        if r['encoding']=='identity' and p.suffix in {'.json','.jsonl','.csv','.log'} and p.stat().st_size>100_000:
            q=Path(str(p)+'.gz')
            with p.open('rb') as i,q.open('wb') as o,gzip.GzipFile(filename='',mode='wb',fileobj=o,mtime=0,compresslevel=3) as g:shutil.copyfileobj(i,g)
            p.unlink();r.update(path=str(q.relative_to(e.root)),bytes=q.stat().st_size,sha256=sha(q),encoding='gzip_exact_original_bytes')
        kept.append(r)
    e.files=kept

sim=load_export('simulation')
extra=BASE/'audit/protected_remote_inputs'
sim.tree(extra,'computation_examples/mevo_acoustic_microbubbles','sim')
sim.file(BASE/'audit/protected_inputs_reconciliation.json','sync/20261022/REMOTE_INPUT_RECONCILIATION.json','shared',force=True)
compact(sim);sim.finish()

exp=load_export('experiment')
v3=Path('/home/lzy/projects/computation_examples/_archives/mevo_vision_f10w_20261017/external_connect/vision_tracking_research/collective_motion_v3')
# Apply the same export policy, plus exclude upstream packaged checkpoints/data.
start=len(exp.files)
exp.tree(v3,'external_connect/vision_tracking_research/collective_motion_v3','vision')
keep=[]
for r in exp.files:
    s=r['path'];p=exp.root/s
    if '/collective_motion_v3/third_party/' in s and (p.suffix in {'.pt','.mat'} or p.name=='fullvideo_delivery.tar.gz') and r['bytes']>5_000_000:
        exp.excluded.append(dict(source=r['source'],path=s,source_bytes=r['source_bytes'],sha256=r['source_sha256'],reason='downloadable_upstream_model_or_demo_dataset'));p.unlink()
    else:keep.append(r)
exp.files=keep
compact(exp);exp.finish()
for e in [sim,exp]:
    shutil.copy2(__file__,e.root/'sync/20261022/finalize_snapshots.py')
    policy=dict(snapshot_date='2026-10-22',branch_type='independent_curated_snapshot_no_imported_development_git_history',
        source_worktrees='read_only',scientific_results='unchanged; no simulation or training executed',
        compressed_files='gzip is lossless; restore_payload.py restores the exact original SHA-256 before reuse',
        retained=['source/config/tests/dependency specifications','scientific reports/failed cases/logs','final figures','critical input data and compact numerical results','manual labels/review scopes/splits and selected trained weights'],
        excluded=['virtual environments/caches/Git internals/vendor installers','full long raw and rendered videos, except one short real input clip','downloadable upstream model weights','regenerable bulk synthetic arrays','superseded unpromoted SAM decoder weights','bulky LF diagnostics and duplicate legacy review/plot arrays'],
        important_external_artifacts='Full experimental recordings and server raw HF trajectories remain external, not a disaster recovery backup. Paths/hashes are retained in EXCLUDED_FILES and historical REMOTE_RAW_INDEX.',
        fine_acoustic_caches='Recomputable fine caches stay on server; exact paths/hashes and generators retained. Coarse caches and response tables included.',
        credentials='No private keys, tokens, or environments are intentionally included. Controller state evidence is archival; no commands are sent.',
        source_identity='Git commit alone does not identify dirty vision/controller edits: COPIED_FILES provides actual source-byte hashes.',
        inherited_manifests='Historical manifests may list omitted files and absolute execution paths; sync/20261022 is the authoritative inventory of this export.')
    write(e.root/'sync/20261022/EXPORT_POLICY.json',policy)
    (e.root/'.gitignore').write_text('.venv*/\n__pycache__/\n.pytest_cache/\n*.pyc\n.env\n')
print('finalization complete')
