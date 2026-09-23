"""Verify local copies against the server manifest; never rewrite audit arrays."""
from pathlib import Path
import json,hashlib,datetime,sys
ROOT=Path(__file__).resolve().parents[1]

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    data=ROOT/'data/remote_results';manifest=data/'OUTPUT_SHA256.json'
    expected=json.loads(manifest.read_text());missing=[];different=[];total=0
    for name,digest in expected.items():
        p=data/name
        if not p.is_file():missing.append(name);continue
        total+=p.stat().st_size
        if sha(p)!=digest:different.append(name)
    result=dict(status='PASS' if not missing and not different else 'FAIL',
                timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                server_output_manifest_sha256=sha(manifest),expected_file_count=len(expected),
                bytes_verified=total,missing_files=missing,mismatched_files=different,
                scope='LOCAL_DOWNLOAD_SHA256_PARITY_WITH_FINAL_SERVER_OUTPUTS')
    (ROOT/'data/delivery_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));return 0 if result['status']=='PASS' else 1

if __name__=='__main__':sys.exit(main())
