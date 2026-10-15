"""WSL delivery readback; write the audit outside the sealed payload."""
from pathlib import Path
import hashlib,json,csv,re,struct,time
R=Path('/home/lzy/projects/temp_storage/clot_github_sync_20261011/clot_clearance');TASK='lps_variational_softening_v3_integration_and_single_pilot_v1';RUN='20261014T232612Z';O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
m=json.loads((O/'PAYLOAD_MANIFEST.json').read_text());assert sha(O/'PAYLOAD_MANIFEST.json')==(O/'PAYLOAD_MANIFEST.sha256').read_text().split()[0]
for name,z in m['files'].items():
 p=R/name;assert p.is_file() and p.stat().st_size==z['bytes'] and sha(p)==z['sha256'],name
# Live initial transfer may have copied pre-archive files. Remove only originals
# whose exact bytes are already verified recoverable in the sealed archive.
removed=[]
for q in (O/'full/rejected_archives').glob('batch_*.json'):
 a=json.loads(q.read_text());assert sha(q.parent/a['archive'])==a['archive_sha256']
 for z in a['files']:
  p=O/'full/rejected'/z['name']
  if p.exists():
   assert str(p.relative_to(R)) not in m['files'] and sha(p)==z['original_sha256'];p.unlink();removed.append(z['name'])
s=json.loads((O/'FINAL_STATUS.json').read_text());doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());actual=Path('/home/lzy/projects/clot_rupture.md');assert sha(actual)==doc['actual_readback_sha256']==doc['after_sha256']
text=(O/('REVIEW_'+TASK.upper()+'.md')).read_text();assert [int(v) for v in re.findall(r'^### (\d+)\.',text,re.M)]==list(range(1,29))
with (O/'V3_NATIVE_PILOT_HISTORY.csv').open() as f:
 rows=list(csv.DictReader(f));assert len(rows)==s['accepted_steps']+1
 for i,row in enumerate(rows):
  assert int(row['step'])==i and int(row['atom_count'])==3488
  assert int(row['attached_nodes'])+int(row['detached_nodes'])==1927
  if i:assert float(row['time_s'])>float(rows[i-1]['time_s'])
assert abs(float(rows[-1]['time_s'])-s['last_accepted_time_s'])<1e-15
pngs={}
for p in F.glob('*.png'):
 with p.open('rb') as f:b=f.read(24)
 assert b[:8]==b'\x89PNG\r\n\x1a\n';w,h=struct.unpack('>II',b[16:24]);assert w>=3000 and h>=1800;pngs[p.name]=[w,h]
assert len(pngs)==len(list(F.glob('*.svg')))==6
for p in [O/'OPEN_RESULTS.html',F/'OPEN_RESULTS.html']:
 for rel in re.findall(r'(?:href|src)="([^"]+)"',p.read_text()):assert (p.parent/rel).resolve().is_file(),(p,rel)
ai=json.loads((O/'V3_SOFTENING_PILOT_AI_RECORD.json').read_text());assert not ai['formal_training_data_eligible'] and ai['CRITICAL_LOAD']=='UNAVAILABLE' and ai['scientific_label']=='EXPLORATORY_UNCALIBRATED'
record=dict(status='PASS',run_id=RUN,manifest_sha256=sha(O/'PAYLOAD_MANIFEST.json'),verified_files=len(m['files']),bytes=m['total_bytes'],scientific_status=s['status'],accepted_steps=s['accepted_steps'],physical_time_s=s['last_accepted_time_s'],WSL_actual_sha256=sha(actual),HTML_targets_valid=True,report_questions=28,PNG_dimensions=pngs,wall_epoch=time.time())
record['redundant_pre_archive_local_copies_removed']=len(removed)
p=R/('LOCAL_DELIVERY_READBACK_'+RUN+'.json');p.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n');print(json.dumps(record,ensure_ascii=False,indent=2))
