"""Verify and restore selected remote result bundles, never overwrite differences."""
import argparse,json,subprocess,hashlib,tarfile,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('bundle',choices=['fullvideo_delivery.tar.gz','models_complete.tar.gz','logs_complete.tar.gz']);a=p.parse_args();index=json.loads((ROOT/'REMOTE_ARTIFACTS.json').read_text());entry=next(x for x in index['bundles'] if x['name']==a.bundle);archive=ROOT/'evidence'/a.bundle;archive.parent.mkdir(exist_ok=True)
    if not archive.exists():subprocess.run(['scp','-i','/home/lzy/.ssh/vast_step3b_ed25519','-P','15283','root@95.3.33.46:'+entry['path'],str(archive)],check=True)
    assert sha(archive)==entry['sha256'],'ARCHIVE_SHA_MISMATCH';dest=ROOT/'evidence/server_bundle' if a.bundle=='logs_complete.tar.gz' else ROOT
    with tarfile.open(archive,'r:gz') as tf:
        for m in tf.getmembers():
            if m.isdir():continue
            assert m.isfile(),'ONLY_REGULAR_RESULT_FILES_ALLOWED';target=dest/m.name;assert target.resolve().is_relative_to(dest.resolve())
            with tf.extractfile(m) as source:
                blob=source.read();h=hashlib.sha256(blob).hexdigest()
                if target.exists():
                    if sha(target)!=h:print('PRESERVE existing differing local file:',target)
                    continue
                target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(blob)
    print('Restored missing files from verified bundle:',a.bundle)
if __name__=='__main__':main()
