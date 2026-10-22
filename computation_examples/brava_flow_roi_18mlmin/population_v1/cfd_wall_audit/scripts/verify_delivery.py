"""Final artifact consistency and preservation checks; no solver or dynamics call."""
from pathlib import Path
import json,hashlib,re
from PIL import Image
root=Path(__file__).resolve().parents[2];audit=root/'cfd_wall_audit';protocol=root/'fixed_window_protocol'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def countrows(p):
 with open(p,'rb') as f:return sum(b.count(b'\n') for b in iter(lambda:f.read(1<<20),b''))-1
badjson=[];jsoncount=0
for folder in [audit,protocol]:
 for p in folder.rglob('*.json'):
  try:json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)));jsoncount+=1
  except Exception as e:badjson.append([str(p),str(e)])
assert not badjson,badjson
counts={name:countrows(audit/name) for name in ['all_cells.csv','wall_adjacent_cells.csv','wall_sample_metrics.csv','wall_nodes.csv','particle_overlap.csv','particle_exposure.csv']}
expected={'all_cells.csv':692313,'wall_adjacent_cells.csv':283313,'wall_sample_metrics.csv':505370,'wall_nodes.csv':50581,'particle_overlap.csv':692313,'particle_exposure.csv':1565};assert counts==expected
figs=[]
for folder,n in [(audit/'figures',12),(protocol/'figures',2)]:
 pngs=sorted(folder.glob('*.png'));assert len(pngs)==n
 for p in pngs:
  im=Image.open(p);im.verify();im=Image.open(p);assert min(im.size)>1000
  assert p.with_suffix('.pdf').stat().st_size>1000
  figs.append(dict(file=str(p.relative_to(root)),pixels=im.size,dpi=im.info.get('dpi')))
spec=json.loads((root/'data/environment.local.json').read_text());inputs=dict(mesh=spec['mesh'],flow=spec['flow'],**spec['boundaries']);inputhash={k:sha(v) for k,v in inputs.items()};assert inputhash==spec['sha256']
links=[]
for p in [audit/'OPEN_RESULTS.html',protocol/'OPEN_RESULTS.html']:
 for target in re.findall(r'(?:href|src)="([^"]+)"',p.read_text()):assert (p.parent/target).exists(),target;links.append(target)
for name in ['history_immutability.json','unchanged_physics_sources.json']:assert json.loads((protocol/'validation'/name).read_text()).get('unchanged',json.loads((protocol/'validation'/name).read_text()).get('all_unchanged'))
log=(protocol/'logs/complete_suite_final.log').read_text();assert '8 passed' in log
empty=json.loads((protocol/'validation/empty_real_environment_fixedwindow_current/MANIFEST.json').read_text());assert empty['physical_time_s']==11.240815660 and empty['actual_source_count']==0 and empty['current_protocol']=='FIXED_WINDOW_POPULATION_EXPERIMENT_V1' and empty['cfd_physical_accuracy_validated_by_protocol'] is False
result=dict(PASS=True,json_files_parsed=jsoncount,json_errors=badjson,csv_rows=counts,figures=figs,html_local_links=len(links),input_hashes_unchanged=True,history_unchanged=True,legacy_physics_unchanged=True,new_protocol_tests=8,large_production_started=False,CFD_executed=False)
(audit/'validation/FINAL_DELIVERY_QC.json').write_text(json.dumps(result,indent=2)+'\n')
for folder in [audit,protocol]:
 files={str(p.relative_to(folder)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(folder.rglob('*')) if p.is_file() and p.name not in ('DELIVERY_MANIFEST.json','final_delivery_qc.log') and '__pycache__' not in p.parts}
 (folder/'DELIVERY_MANIFEST.json').write_text(json.dumps(dict(files=files,file_count=len(files),data_role='READ_ONLY_AUDIT' if folder==audit else 'FIXED_WINDOW_PROTOCOL_AND_SMALL_VALIDATION',code_version=json.loads((protocol/'CODE_VERSION.json').read_text())),indent=2)+'\n')
print(json.dumps(result,indent=2))
