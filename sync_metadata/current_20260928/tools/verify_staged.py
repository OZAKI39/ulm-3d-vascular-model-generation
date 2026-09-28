"""Compare the exact staged Git tree with the content inventory, without filters."""
from pathlib import Path
import csv,hashlib,json,os,subprocess
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260928'
def blob(p):
    data=os.readlink(p).encode() if p.is_symlink() else p.read_bytes()
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
def main():
    raw=subprocess.check_output(['git','ls-files','--stage','-z'],cwd=D)
    staged={}
    for line in raw.split(b'\0'):
        if not line:continue
        desc,path=line.split(b'\t',1);mode,oid,stage=desc.decode().split();assert stage=='0'
        staged[path.decode()]={'mode':mode,'oid':oid}
    inventory={r['path']:r for r in csv.DictReader((M/'snapshot_manifest.csv').open())}
    failures=[]
    for path,r in inventory.items():
        if path not in staged:failures.append((path,'missing from index'));continue
        if staged[path]['oid']!=r['git_blob_oid'] or staged[path]['mode']!=r['mode']:failures.append((path,'staged content/mode mismatch'))
    extras=sorted(set(staged)-set(inventory))
    for path in extras:
        assert path.startswith('sync_metadata/current_20260928/'),path
        if blob(D/path)!=staged[path]['oid']:failures.append((path,'generated file differs from index'))
    result={'PASS':not failures,'staged_files_and_links':len(staged),'manifest_entries_verified':len(inventory),
        'generated_extras_verified':extras,'failures':failures,'scientific_filters_disabled':True}
    (M/'validation/staged_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2));assert result['PASS']
if __name__=='__main__':main()
