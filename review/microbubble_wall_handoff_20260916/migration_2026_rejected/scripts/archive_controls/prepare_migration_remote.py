from pathlib import Path
import json,hashlib,shutil,difflib,subprocess
R=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244');W=Path('/workspace/lammps_migration/new_20260916_091244')
OLD={'engine':'lammps_particle_engine_20260915_214759','coupling':'palabos_lammps_coupling_v0_20260915_224019','passive':'passive_transport_v0_20260915_233516'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
provenance={};patches=[]
for stage,name in OLD.items():
 src=Path('/workspace/microbubble_lammps/results')/name;dst=R/'migration'/stage;assert not dst.exists();dst.mkdir(parents=True)
 manifest=src/'FROZEN_INPUT_SHA256SUMS'
 for line in manifest.read_text().splitlines():
  h,n=line.split('  ',1);assert sha(src/n)==h,n
  f=dst/n;f.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src/n,f)
 shutil.copy2(manifest,dst/'BASELINE_FROZEN_INPUT_SHA256SUMS')
 for sub in ['src','contracts']:
  shutil.copytree(src/sub,dst/sub,dirs_exist_ok=True)
 if stage=='coupling':
  for n in ['COORDINATE_NODE_CHECK.csv']:
   shutil.copy2(src/'validation'/n,dst/'validation'/n)
 if stage in ['coupling','passive']:
  n='src/coupled_lammps_main.cpp' if stage=='coupling' else 'src/passive_main.cpp'
  f=dst/n;before=f.read_text();after=before.replace('(*lmp->modify->fix_map)["frozen/flow/drag"]=&make_drag','Modify::fix_styles().set_plugin("frozen/flow/drag", &make_drag)').replace('(*l->modify->fix_map)["sonovue/passive/flow"]=&factory','Modify::fix_styles().set_plugin("sonovue/passive/flow", &factory)')
  assert before!=after
  f.write_text(after);patches.append({'path':str(f),'old_sha256':hashlib.sha256(before.encode()).hexdigest(),'new_sha256':sha(f),'reason':'2Sep2026 global CreatorRegistry API replaces instance fix_map; only factory registration changes','diff':''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='old/'+n,tofile='new/'+n))})
  cm=dst/'CMakeLists.txt';text=cm.read_text();assert text== (src/'CMakeLists.txt').read_text()
 # All baseline numerical fields, states, populations and gates must remain byte-identical.
 mismatches=[]
 for line in manifest.read_text().splitlines():
  h,n=line.split('  ',1)
  if sha(dst/n)!=h and not n.endswith(('passive_main.cpp','coupled_lammps_main.cpp')):mismatches.append(n)
 assert not mismatches
 provenance[stage]={'historical_root':str(src),'historical_manifest_sha256':sha(src/'SHA256SUMS'),'copied_input_manifest_sha256':sha(manifest),'math_and_inputs_unchanged':True}
 lines=[sha(f)+'  '+str(f.relative_to(dst)) for f in sorted(dst.rglob('*')) if f.is_file() and f.name!='FROZEN_INPUT_SHA256SUMS']
 (dst/'FROZEN_INPUT_SHA256SUMS').write_text('\n'.join(lines)+'\n')
(R/'provenance/MIGRATION_INPUT_PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
(R/'provenance/API_COMPATIBILITY_PATCH.json').write_text(json.dumps(patches,indent=2)+'\n')
S=W/'upstream/lammps';old=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759/upstream/lammps');diffs=[]
for name in ['pair_lubricate_poly.cpp','pair_lubricateU_poly.cpp']:
 a=(old/'src/COLLOID'/name).read_text();b=(S/'src/COLLOID'/name).read_text();diffs.append(''.join(difflib.unified_diff(a.splitlines(True),b.splitlines(True),fromfile='22Jul2025_update6/'+name,tofile='2Sep2026/'+name)))
(R/'provenance/upstream_lubrication/OLD_TO_TARGET.diff').write_text('\n'.join(diffs))
for n in ['modify.h','creator_registry.h']:
 shutil.copy2(S/'src'/n,R/'provenance'/n)
print(json.dumps({'status':'PREPARED','stages':list(OLD),'compatibility_patches':patches},indent=2))
