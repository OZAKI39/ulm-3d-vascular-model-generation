#!/usr/bin/env python3
"""Record every retained remote output and CSV point summary for review/sync."""
from pathlib import Path
import argparse,csv,gzip,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_diagnostics import basin_rows
from particle_3d.particle82_provenance import atomic_json,sha256,require_remote

p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--points',required=True);a=p.parse_args()
r=Path(a.root);points=Path(a.points);h=json.loads((r/'host_provenance.json').read_text());require_remote(h,h['hostname'])
(r/'data').mkdir(exist_ok=True)
for label,folder in [('point_baseline',points),('point_refinement',r/'point_refinement_100000'),('original_proposal_points',r/'admission/point_paths')]:
 rows=basin_rows(folder);columns=list(dict.fromkeys(k for row in rows for k in row))
 with gzip.open(r/'data'/(label+'_basins.csv.gz'),'wt',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)
files=[]
for label,folder in [('full_audit',r),('point_basin_100000',points)]:
 for f in sorted(folder.rglob('*')):
  if not f.is_file() or f.is_relative_to(r/'repo') or f.name in ['REMOTE_RAW_INVENTORY.json','REMOTE_RAW_SHA256SUMS.txt']:continue
  files.append(dict(collection=label,relative_path=str(f.relative_to(folder)),remote_path=str(f),bytes=f.stat().st_size,sha256=sha256(f)))
atomic_json(r/'REMOTE_RAW_INVENTORY.json',dict(hostname=h['hostname'],endpoint=h['ssh_endpoint'],files=files,
 total_bytes=sum(f['bytes'] for f in files),status='PRESERVED_AT_RECORDED_REMOTE_PATHS; LOCAL_COPY_STATUS_REQUIRES_LOCAL_SYNC_INVENTORY'))
(r/'REMOTE_RAW_SHA256SUMS.txt').write_text(''.join(f'{f["sha256"]}  {f["remote_path"]}\n' for f in files))
print(len(files),sum(f['bytes'] for f in files))
