"""Check exact Git blobs and record a snapshot manifest before publication."""
from pathlib import Path
import collections,csv,hashlib,json,os,subprocess
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260929'
assert json.loads((M/'validation/snapshot_validation.json').read_text())['PASS']
index={}
for line in subprocess.check_output(['git','ls-files','--stage','-z'],cwd=D).split(b'\0'):
    if not line:continue
    desc,path=line.split(b'\t',1);mode,oid,stage=desc.decode().split();assert stage=='0'
    index[path.decode()]=(mode,oid)
manifest=[];failures=[]
for rel,(mode,oid) in index.items():
    assert mode!='160000','Unexpected submodule '+rel
    p=D/rel
    if mode=='120000':
        data=os.readlink(p).encode();size=len(data)
        digest=hashlib.sha256(data).hexdigest();blob=hashlib.sha1(b'blob '+str(size).encode()+b'\0'+data).hexdigest()
    else:
        size=p.stat().st_size;s256=hashlib.sha256();s1=hashlib.sha1(b'blob '+str(size).encode()+b'\0')
        with p.open('rb') as f:
            while block:=f.read(1024*1024):s256.update(block);s1.update(block)
        digest=s256.hexdigest();blob=s1.hexdigest()
    if blob!=oid:failures.append(rel)
    assert size<100*1024**2,rel
    manifest.append(dict(path=rel,mode=mode,bytes=size,sha256=digest,git_blob=oid))
for row in csv.DictReader((M/'source_files.csv').open()):assert row['destination'] in index,row['destination']
for row in json.loads((M/'server_file_map.json').read_text()):
    if 'destination' in row:assert row['destination'] in index,row['destination']
with (M/'snapshot_manifest.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['path','mode','bytes','sha256','git_blob']);writer.writeheader();writer.writerows(manifest)
result=dict(PASS=not failures,indexed_files_and_links=len(index),logical_bytes=sum(r['bytes'] for r in manifest),
    maximum_blob_bytes=max(r['bytes'] for r in manifest),Git_blob_mismatches=failures,
    root_counts=dict(collections.Counter(Path(p).parts[0] for p in index)),
    manifest_scope='All indexed files before this manifest and its verification result are added; these two files do not self-reference')
(M/'validation/index_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));assert result['PASS']
