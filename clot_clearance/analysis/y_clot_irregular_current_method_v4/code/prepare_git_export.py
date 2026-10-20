"""Local archive selection only; no Git mutation and no network push."""
from pathlib import Path
import json,gzip,hashlib,datetime
R=Path(__file__).resolve().parents[3];TASK='y_clot_irregular_current_method_v4';RUN='20261020T041100Z';D=R/'analysis'/TASK;O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 assert json.loads((D/'FINAL_STATUS.json').read_text())['delivery_validation']=='PASS'
 records=[];selected={};omitted={}
 excluded={'GIT_EXPORT_MANIFEST.json','GIT_STAGE_PATHS.txt','GIT_COMPRESSED_EVIDENCE.json'}
 for base in [D,O,F]:
  for p in sorted(base.rglob('*')):
   if not p.is_file() or '__pycache__' in p.parts or p.name in excluded or p.name.endswith('.rendering.mp4'):continue
   key=str(p.relative_to(R));size=p.stat().st_size
   if p.name=='LOAD_RENDER_PREVIEW.png':
    omitted[key]=dict(bytes=size,reason='Obsolete early-load rendering preview; final actual load field delivered');continue
   if size>=32*2**20:
    omitted[key]=dict(bytes=size,sha256=sha(p),reason='Large raw artifact retained locally and on server');continue
   if p.is_relative_to(O) and p.suffix in ['.json','.csv','.log','.txt','.data'] and size>2**20:
    target=p.with_name(p.name+'.gz');raw=p.read_bytes();packed=gzip.compress(raw,mtime=0);assert gzip.decompress(packed)==raw
    target.write_bytes(packed);records.append(dict(source=str(p.relative_to(O)),source_bytes=len(raw),source_sha256=hashlib.sha256(raw).hexdigest(),compressed=str(target.relative_to(O)),compressed_bytes=len(packed),compressed_sha256=sha(target)))
    omitted[key]=dict(bytes=size,sha256=sha(p),reason='Losslessly gzip compressed; source file remains on local disk')
    selected[str(target.relative_to(R))]=dict(bytes=target.stat().st_size,sha256=sha(target))
   else:selected[key]=dict(bytes=size,sha256=sha(p))
 compression=O/'GIT_COMPRESSED_EVIDENCE.json';compression.write_text(json.dumps(dict(records=records),indent=2)+'\n')
 selected[str(compression.relative_to(R))]=dict(bytes=compression.stat().st_size,sha256=sha(compression))
 manifest=D/'GIT_EXPORT_MANIFEST.json';manifest.write_text(json.dumps(dict(status='PREPARED_AND_HASHED',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),selected=selected,omitted=omitted,restore=f'python analysis/{TASK}/code/restore_export.py',no_github_push=True),indent=2)+'\n')
 paths=sorted([*selected,str(manifest.relative_to(R))]);(D/'GIT_STAGE_PATHS.txt').write_bytes(b'\0'.join(p.encode() for p in paths)+b'\0')
 print('PREPARED',len(selected),'files;',len(records),'lossless text compressions;',len(omitted),'omissions')
if __name__=='__main__':main()
