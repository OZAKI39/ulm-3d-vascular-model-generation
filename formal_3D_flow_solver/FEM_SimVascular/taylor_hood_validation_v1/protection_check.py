"""Final read-only hash/git comparison. Never attempts automatic restoration."""
from common import *
import subprocess,os,time
from concurrent.futures import ThreadPoolExecutor

def compare(item):
    path,before=item;p=Path(path)
    if not p.is_file():return dict(path=path,before=before,after=None,reason='missing or no longer a regular file')
    after=dict(size=p.stat().st_size,sha256=sha(p))
    return dict(path=path,before=before,after=after) if any(before[k]!=after[k] for k in ('size','sha256')) else None

def main():
    start=time.time();manifest=json.loads((REPORT/'baseline_protection_manifest.json').read_text())
    with ThreadPoolExecutor(max_workers=6) as pool:changes=[v for v in pool.map(compare,manifest['files'].items()) if v]
    before=json.loads((REPORT/'logs/git_before.json').read_text());after={};git_changes=[]
    for repo,commands in before.items():
        after[repo]={}
        for command in commands:
            p=subprocess.run(['git','-C',repo,*command.split()],capture_output=True,text=True,env=dict(os.environ,GIT_OPTIONAL_LOCKS='0'))
            row=dict(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr);after[repo][command]=row
            if row['stdout']!=commands[command]['stdout']:
                git_changes.append(dict(repository=repo,command=command,before=commands[command]['stdout'],after=row['stdout']))
    dump(REPORT/'logs/git_after.json',after)
    result=dict(status='PASS' if not changes else 'PROTECTED_FILE_MODIFICATION_DETECTED',
        protected_files_checked=len(manifest['files']),protected_files_changed=len(changes),changed_files=changes,
        git_observed_differences=git_changes,git_difference_note='New task directories may appear as untracked; existing HEAD, index and files are separately SHA-protected.',
        elapsed_seconds=time.time()-start,restoration_attempted=False)
    dump(REPORT/'data/local_protection_final.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['git_observed_differences','changed_files']}),flush=True)

if __name__=='__main__':main()
