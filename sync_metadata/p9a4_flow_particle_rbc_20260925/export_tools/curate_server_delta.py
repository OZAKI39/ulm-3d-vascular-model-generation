from pathlib import Path
import json,gzip,hashlib,subprocess,collections
A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination']);M=D/'sync_metadata'/C['snapshot']
x=json.loads((A/'audit/server_delta_inventory.json').read_text());tracked=set(subprocess.check_output(['git','ls-files','-z'],cwd=D).decode().split('\0'))
def vendor(p):
 parts=Path(p).parts
 return bool(set(parts)&{'petsc-3.19.6','petsc_source','palabos','externalLibraries','tbb'}) or '/toolchain_env/installer/' in p
removed=[]
for r in x['files']:
 target=r.get('representation',{}).get('path','') if isinstance(r.get('representation'),dict) else r['target']
 if vendor(r['source']) or vendor(target):
  if r['status']=='EXPORTED' and r['target'] not in tracked:
   p=D/r['target']
   if p.exists():p.unlink();removed.append(r['target'])
  r['status']='EXCLUDED_FROM_DELTA';r['reason']='UNMODIFIED_VENDOR_SOURCE_OR_TOOLCHAIN_INSTALLER';r.pop('representation',None)
x['exported_archive_scope']='Transport archive contains extra vendor files pruned after inspection; this manifest describes final repository selection.'
x['removed_vendor_files']=len(removed)
# Check every content reference against actual published bytes, with gzip handled explicitly.
cache={};invalid=[]
for r in x['files']:
 if r['status']!='REPRESENTED_BY_EXISTING_ARTIFACT':continue
 rep=r['representation'];p=D/rep['path'];key=(rep['path'],rep['representation'])
 if key not in cache:
  if not p.is_file():cache[key]=None
  else:
   opener=gzip.open if rep['representation']=='GZIP_OF_IDENTICAL_BYTES' else open
   with opener(p,'rb') as f:cache[key]=hashlib.file_digest(f,'sha256').hexdigest()
 if cache[key]!=r['sha256']:invalid.append(dict(source=r['source'],target=rep['path'],expected=r['sha256'],actual=cache[key]))
(A/'audit/server_reference_validation.json').write_text(json.dumps(dict(unique_targets=len(cache),invalid=invalid,removed_vendor_files=len(removed)),indent=2)+'\n')
if invalid:raise RuntimeError('Invalid server content references: '+str(invalid[:4]))
p=M/'server_delta_inventory.json.gz'
with p.open('xb') as raw:
 with gzip.GzipFile(fileobj=raw,mode='wb',filename='',mtime=0,compresslevel=6) as f:f.write((json.dumps(x,indent=2)+'\n').encode())
(M/'server_delta_inventory.json').unlink() # new export copy only; full transport record remains in local archive
summary=dict(counts=dict(collections.Counter(r['status'] for r in x['files'])),new_files=sum(r['status']=='EXPORTED' for r in x['files']),new_bytes=sum(r['size'] for r in x['files'] if r['status']=='EXPORTED'),vendor_files_excluded=len(removed),content_references_verified=True,inventory_compressed_bytes=p.stat().st_size)
(M/'server_delta_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
