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
shard_records=[]
for f in sorted(r.glob('**/shards/shard_*.json')):
 if f.is_relative_to(r/'repo'):continue
 row=json.loads(f.read_text());folder=f.parent.parent
 if 'start_time' not in row:
  receipts=[json.loads((folder/p).read_text()) for p in row['individual_receipts']]
  row['start_time']=min(x['start_time'] for x in receipts)
  row['start_time_source']='DERIVED_FROM_ORIGINAL_CONTEMPORANEOUS_PER_ID_RECEIPTS; ORIGINAL_INDEX_UNMODIFIED'
 else:row['start_time_source']='ORIGINAL_SHARD_EXECUTION_RECORD'
 row.update(original_shard_path=str(f),original_shard_sha256=sha256(f))
 assert row['start_time']<=row['end_time'] and row['hostname']==h['hostname']
 shard_records.append(row)
atomic_json(r/'data/SHARD_PROVENANCE_COMPLETENESS.json',dict(all_pass=True,shards=shard_records,
 note='All primary trajectory receipts already contain start/end. Early benchmark secondary indices are supplemented transparently from those original receipts.'))
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
