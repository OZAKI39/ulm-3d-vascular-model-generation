from pathlib import Path
import json,hashlib,subprocess,shutil
R=Path(__file__).resolve().parents[1];W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work']);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();res=[]
for b in json.loads((R/'provenance/BASELINES_BEFORE.json').read_text()):
 root=Path(b['root']);p=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=root,capture_output=True,text=True);name=root.parent.name+'_'+root.name;(R/'logs'/f'IMMUTABLE_{name}.log').write_text(p.stdout+p.stderr);res.append(dict(root=str(root),manifest_unchanged=sha(root/'SHA256SUMS')==b['manifest_sha256'],all_manifest_files_pass=p.returncode==0));assert res[-1]['manifest_unchanged'] and res[-1]['all_manifest_files_pass'],root
provenance=json.loads((R/'OPEN_SOURCE_TOOL_PROVENANCE.json').read_text())
for repo in provenance:
 expected=json.loads((R/'provenance'/repo['source_sha_manifest']).read_text());root=Path(repo['local_source']);actual={p:sha(root/p) for p in expected};assert actual==expected,repo['name'];(R/'provenance'/f"{repo['name']}_SOURCE_AFTER.json").write_text(json.dumps(actual,indent=2)+'\n');res.append(dict(repo=repo['name'],tracked_source_files=len(expected),status='PASS',source_modified=False));repo['use']='PRIMARY_CAPABILITY_PASS' if repo['name']=='Pecnut' else 'CAPABILITY_AUDIT_ONLY'
 if repo['name']=='RMBW':repo['backend']='multiple native backends available; no runtime backend invoked in this audit'
 if repo['name']=='RigidBodyIB':repo['backend']='standalone C/C++ marker mobility; capability only'
# Native mathematical files in the execution copy are byte-identical to upstream.
orig=W/'upstream/stokesian-dynamics/stokesian_dynamics';copy=W/'runtime/stokesian_dynamics';native={str(p.relative_to(orig)):sha(p) for p in orig.rglob('*.py')};assert all(sha(copy/p)==h for p,h in native.items());res.append(dict(runtime_native_Python_files=len(native),status='PASS'))
for name in ['rigid_math.cpp','rigid_math.hpp','frozen_flow.hpp']:
 original=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['pair_baseline'])/'src'/name;assert sha(original)==sha(R/'src/pair_baseline_copy'/name)
assert sha(R/'contracts/FIXED_MULTIBLOB_WALL_FEASIBILITY_CONTRACT.json')=='781deaa0af68f84196f0dda01db078ac30f64a9aafed065bc9bbe89483702a58'
assert sha(R/'input/RMBW_WALL_RESISTANCE_TABLE_V0.h5')=='71e62ef6bc22edc66ae0c4d1dcc224378c2f6c5d3848578f3fbedd2cd9a43f98'
assert sha(R/'input/closed_geometry_m.stl')=='840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb'
(R/'OPEN_SOURCE_TOOL_PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n');(R/'validation/IMMUTABILITY_AFTER.json').write_text(json.dumps(dict(status='PASS',results=res),indent=2)+'\n');print('IMMUTABILITY_PASS',len(res))
