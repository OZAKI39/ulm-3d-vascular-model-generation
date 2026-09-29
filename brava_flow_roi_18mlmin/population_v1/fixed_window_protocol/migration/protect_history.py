from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parents[2]
paths=[root/'runs/realization_0000_contactfix',root/'analysis',root/'candidate_0',root/'candidate_15',root/'reports']
paths += [Path('/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/src/particle_3d')]
exclude={root/'fixed_window_protocol'}
result={}
for base in paths:
 for p in sorted(base.rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts:continue
  if base.name=='particle_3d' and ('population_v1' in p.parts or p.suffix not in ('.py','.cu','.cpp')):continue
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1<<20),b''):h.update(b)
  result[str(p)]=h.hexdigest()
output=root/'fixed_window_protocol/validation'/sys.argv[1]
output.write_text(json.dumps(result,indent=2)+'\n');print(len(result),output)
