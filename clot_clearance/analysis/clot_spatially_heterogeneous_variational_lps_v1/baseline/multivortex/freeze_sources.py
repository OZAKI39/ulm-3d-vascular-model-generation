from common import *
import subprocess,difflib
assert json.loads((O/'KEY_CHECKS.json').read_text())['status']=='PASS'
assert not (O/'SOURCE_FREEZE.json').exists()
surface=json.loads((O/'CURRENT_EXPOSED_SURFACE_SPEC.json').read_text())
surface['current_occlusion']='Gap-cell centroid, portal centroid and two outward face probes at 1/3 um must lie outside every current fully supported tetra (KD-tree broad phase and exact barycentric point containment). Occluded candidates remain unresolved.'
surface['geometric_limitations'][-1]='Finite point-probe visibility cannot certify a continuous liquid corridor or resolve subgrid throats/self-contact. No infiltration PDE. Report this approximation and probe rejections.'
dump(O/'CURRENT_EXPOSED_SURFACE_SPEC.json',surface)
paths=[p for p in A.glob('*.py') if p.name not in ['freeze_sources.py']]
dependencies=[R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',PA/'native.py',LIB]
inputs=[p for p in (O/'provenance').iterdir() if p.is_file() and p.name!='new_surface_and_controller.diff']+[O/q for q in ['MATERIAL_AND_SPECIMEN_FREEZE.json','CURRENT_EXPOSED_SURFACE_SPEC.json','WORK_MATCH_SPEC.json','SHEAR_PARAMETERIZATION.json','MULTIFRAGMENT_OBJECTIVE.json','KEY_CHECKS.json']]
dump(O/'SOURCE_FREEZE.json',dict(git_sha=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip(),implementation_sha256={str(p.relative_to(R)):sha(p) for p in paths+dependencies},input_sha256={str(p.relative_to(O)):sha(p) for p in inputs},library_sha256=sha(LIB),main_trial_limit=2,main_trial_names=CASES,history='Read-only accepted prior files; all new code isolated',surface_scope='Conservative supported-simplex / dilated fluid-gap proxy; no continuum completeness claim'))
diff=[]
for p in paths:
    old=OLD_A/p.name
    diff.extend(difflib.unified_diff(old.read_text().splitlines(True) if old.exists() else [],p.read_text().splitlines(True),fromfile=str(old),tofile=str(p)))
(O/'provenance/new_surface_and_controller.diff').write_text(''.join(diff))
print('FROZEN',sha(O/'SOURCE_FREEZE.json'))
