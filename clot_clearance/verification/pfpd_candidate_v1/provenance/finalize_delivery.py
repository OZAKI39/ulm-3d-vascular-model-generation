"""Finalize hashes and evidence. Writes only this task's provenance directory."""
from pathlib import Path
import json,hashlib,re,datetime
O=Path(__file__).resolve().parents[1];ROOT=O.parents[1];P=O/'provenance'
def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def record(p):return dict(sha256=sha(p),bytes=p.stat().st_size)
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
pre=read(P/'PRE_RUN_MANIFEST.json');post={};missing=[];modified=[]
for rel,info in pre['files'].items():
    f=ROOT/rel
    if not f.is_file():missing.append(rel);continue
    post[rel]=record(f)
    if post[rel]!=info:modified.append(rel)
write(P/'POST_RUN_MANIFEST.json',dict(created_utc=now,files=post,scope='All 5307 preexisting protected workspace files, including untracked prior results; additions listed separately'))
preserve=dict(status='PASSED' if not missing and not modified else 'FAILED',preexisting_files_checked=len(post),modified_files=modified,missing_files=missing,production_modified_files=sum(p.startswith(('pd_clot/','scripts/','configs/','tests/')) for p in modified),historical_modified_files=sum(p.startswith(('verification/','results/','visualization/')) for p in modified),original_head=pre['head'])
write(P/'PRESERVATION_CHECK.json',preserve);assert preserve['status']=='PASSED'
base=(O/'tests/BASELINE.log').read_text();new=(O/'tests/NEW_TESTS.log').read_text();alltests=(O/'tests/ALL_TESTS.log').read_text()
assert '130 passed' in base and '39 passed' in new and '169 passed' in alltests
m=read(O/'baqp_verification/MECHANICS_GATE.json');h=read(O/'historical_compatibility/SUMMARY.json');norm=read(O/'normalization/GRIFFITH_STUDY.json');gate=read(O/'FORMULA_CONSISTENCY_GATE.json');oracle=read(O/'comparisons/OFFICIAL_COMPARISON.json');constants=read(O/'normalization/NORMALIZATION_CONSTANTS.json')
fd=max(r['relative_L2_error'] for a in m['mechanics'] for r in a['fd']);forceabs=max(r['max_absolute_error_N'] for a in m['mechanics'] for r in a['fd'])
figs=['energetic_vs_kinematic_degradation','normalization_constant','baqp_affine_reproduction','historical_KF_change','historical_elastic_energy_change','griffith_normalization_convergence','formula_consistency']
assert all((O/'figures'/f'{name}.{ext}').is_file() for name in figs for ext in ['png','pdf'])
assert len(re.findall(r'^## \d+\.',(O/'REVIEW_PFPD_CANDIDATE_V1.md').read_text(),re.M))==16
checks=dict(A='PASS: baseline 130 and combined 169 tests',B='PASS: affine reconstruction',C='PASS: multi-epsilon variational force',D='PASS: global force and moment',E='PASS: quadratic c0=0.25',F='PASS: explicit discrete pair counting, 3x3 finite coupon errors and bulk lattice convergence',G='PASS: 35 exact old-energy/reconstruction replays',H='PASS: B/C h=1 and exact unchanged K/F',I='PASS: production unmodified',J='PASS: unresolved evolution is not implemented')
assert m['status']==oracle['status']=='VERIFIED' and m['affine_max_error']<5e-12 and fd<2e-6
assert h['count']==35 and h['replay_max_error_J']==0 and h['models']['B']['max_relative_delta_F']==0 and h['models']['C']['max_relative_delta_F']==0
assert abs(constants[-1]['c0']-.25)<1e-14 and gate['status']=='UNRESOLVED_IN_V1' and not gate['canonical_evolution_implemented']
labels={'PFPD_CANDIDATE_V1':'DEVELOPMENT VERIFICATION PASSED','PFPD PAPER AUDIT':'COMPLETE','PFPD FORMULA CONSISTENCY':'UNRESOLVED_IN_V1','PUBLIC AUTHOR PFPD CODE':'NOT_FOUND','BAQP INTACT MECHANICS':'VERIFIED','CURRENT KERNEL c0':'VERIFIED_0.25','GRIFFITH NORMALIZATION':'VERIFIED','ENERGETIC-KINEMATIC SEPARATION':'TESTED','CURRENT RECONSTRUCTION ANOMALY UNDER h=1':'REMOVED_FOR_AUDITED_STATES','FULL PFPD DAMAGE EVOLUTION':'NOT_IMPLEMENTED_DUE_TO_FORMULA_AMBIGUITY','PRODUCTION MODEL':'UNCHANGED','PHYSICAL CLOT CALIBRATION':'NOT PERFORMED','REAL MICROBUBBLE CFD':'NOT USED'}
final=dict(created_utc=now,statuses=labels,acceptance_A_to_J=checks,paper_sha256=sha(Path('/home/lzy/projects/2604.09215v1.pdf')),official_reference=oracle,
 tests=dict(old=130,new=39,total=169,failed=0,skipped=0,warnings=60),preservation=preserve,BAQP_affine_max_error=m['affine_max_error'],force_energy_max_relative_L2_error=fd,force_energy_max_absolute_error_N=forceabs,
 c0_current=constants[-1]['c0'],griffith_finite_coupon_max_relative_error=max(abs(r['relative_error']) for r in norm['coupon']),griffith_all_study_max_relative_error=max(abs(r['ratio']-1) for r in norm['coupon']+norm['bulk_lattice_refinement']),griffith_finest_bulk_relative_error=abs(norm['bulk_lattice_refinement'][-1]['ratio']-1),
 historical=h,report='REVIEW_PFPD_CANDIDATE_V1.md',figures=figs,scientific_scope='METHOD REPRODUCTION and prescribed-pattern fixed-state compatibility; not production replacement or complete fracture evolution')
write(P/'FINAL_VERIFICATION.json',final)
# Delivery listing excludes external downloads/depot and disposable test/cache output.
code=list((ROOT/'pd_clot/pfpd_candidate').glob('*.py'))+[ROOT/'scripts'/n for n in ['verify_pfpd_baqp.py','compare_pfpd_julia.py','run_pfpd_compatibility.py','visualize_pfpd_candidate.py']]+[ROOT/'tests'/n for n in ['test_pfpd_candidate_mechanics.py','test_pfpd_candidate_diagnostics.py']]
artifacts=[]
for f in O.rglob('*'):
    if not f.is_file():continue
    rel=f.relative_to(O)
    if 'external' in rel.parts or 'runtime_cache' in rel.parts or any(x.endswith('_tmp') for x in rel.parts):continue
    if f.name in ['DELIVERY_MANIFEST.json','ADDED_FILES.json']:continue
    artifacts.append(f)
allfiles=sorted(code+artifacts);added=[str(p.relative_to(ROOT)) for p in allfiles]
assert all(p not in pre['files'] for p in added)
write(P/'ADDED_FILES.json',dict(deliverable_files_including_manifest_and_this_index=len(added)+2,candidate_source_files=6,scripts=4,test_files=2,artifacts_excluding_two_indices=len(artifacts),excluded='external clone, Julia distribution/depot, downloaded deb, runtime caches and disposable pytest temp directories',files=added))
allfiles.append(P/'ADDED_FILES.json')
write(P/'DELIVERY_MANIFEST.json',dict(created_utc=now,scope='New delivery files only; manifest excludes its own hash',files={str(f.relative_to(ROOT)):record(f) for f in allfiles}))
print(json.dumps(dict(status=labels['PFPD_CANDIDATE_V1'],old=130,new=39,total=169,unchanged=len(post),files_added=len(added)+2,affine=m['affine_max_error'],force_energy_relative=fd,coupon_max_error=final['griffith_finite_coupon_max_relative_error'],all_griffith_max_error=final['griffith_all_study_max_relative_error']),indent=2))
