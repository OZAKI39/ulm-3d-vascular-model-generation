"""Inventory the completed audit only; never start the solver."""
from common import *
exclude={'FILE_MANIFEST.json','SHA256SUMS','DELIVERY_VERIFICATION.json','LOCAL_DELIVERY_SHA256.json'}
records=[]
for folder in [O,F]:
 for p in sorted(folder.rglob('*')):
  if p.is_file() and p.name not in exclude:
   records.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)})
summary={}
for label,folder in [('verification',O),('results',F)]:
 prefix=str(folder.relative_to(R))+'/'
 subset=[r for r in records if r['path'].startswith(prefix)];summary[label]={'files':len(subset),'bytes':sum(r['bytes'] for r in subset)}
m={'run_id':C['run'],'path_base':'Project root: server /workspace/y_clot/project or local clot_clearance','files':len(records),'bytes':sum(r['bytes'] for r in records),'by_directory':summary,'excluded_metadata':sorted(exclude),'note':'Manifest/checksums/transfer-verification metadata excluded to avoid self-reference; includes formal and failed setup evidence, reproduction code, document before/after/diff and all figures','records':records}
dump(O/'FILE_MANIFEST.json',m)
lines=[r['sha256']+'  '+r['path'] for r in records];lines.append(sha(O/'FILE_MANIFEST.json')+'  '+str((O/'FILE_MANIFEST.json').relative_to(R)));(O/'SHA256SUMS').write_text('\n'.join(lines)+'\n')
print(json.dumps({'files':m['files'],'bytes':m['bytes'],'by_directory':summary}))
