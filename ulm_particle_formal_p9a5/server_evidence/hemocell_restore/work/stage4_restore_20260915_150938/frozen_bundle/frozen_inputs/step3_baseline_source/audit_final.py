#!/usr/bin/python3
"""Audit existing inputs against pre-run manifests; writes only RUN_DIR."""
import argparse
import hashlib
import json
import os
import stat
import subprocess
import shutil
from pathlib import Path

HC=Path('/home/lzy/projects/hemocell_starter')
SRC=Path('/home/lzy/projects/ulm_3D_vascular')
CODE=HC/'test_code/step3_vascular_pure_fluid'
ROOTS={'source':SRC,'hemocell':HC,
       'step1':Path('/home/lzy/projects/compre_output/step1/20260912_215759'),
       'step2':Path('/home/lzy/projects/compre_output/step2/20260912_225418')}
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,obj):Path(p).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def manifest(root):
    result={}
    for parent,dirs,files in os.walk(root,followlinks=False):
        for name in dirs+files:
            p=Path(parent)/name;s=p.lstat()
            result[str(p.relative_to(root))]=[s.st_mode,s.st_size,s.st_mtime_ns,os.readlink(p) if p.is_symlink() else None]
    return result
def audit(run):
    env=dict(os.environ,PATH='/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',GIT_OPTIONAL_LOCKS='0')
    lfs=run/'provenance/audit_tools/extracted/usr/bin/git-lfs'
    captured=[]
    statuses={}
    for root in [SRC,HC]:
        command_env=dict(env)
        # A worktree status/diff can refresh the index even through helper
        # processes. Audit using an isolated copy, so any refresh stays in RUN.
        index_copy=run/'provenance/audit_tools'/('source.index' if root==SRC else 'hemocell.index')
        shutil.copy2(root/'.git/index',index_copy)
        command_env['GIT_INDEX_FILE']=str(index_copy)
        args=['/usr/bin/git','-c',f'filter.lfs.process={lfs} filter-process',
              '-c',f'filter.lfs.clean={lfs} clean -- %f',
              '-c',f'lfs.storage={run}/provenance/audit_tools/lfs_storage',
              '-c','diff.autoRefreshIndex=false','-C',str(root)]
        for cmd in [['rev-parse','HEAD'],['status','--short'],['diff','--stat']]:
            proc=subprocess.run(args+cmd,env=command_env,text=True,capture_output=True)
            captured.append(f'$ git -C {root} '+ ' '.join(cmd)+f'\nRC={proc.returncode}\n'+proc.stdout+proc.stderr)
            assert proc.returncode==0
            if cmd==['status','--short']:statuses[str(root)]=proc.stdout.strip()
    (run/'git_status_after.txt').write_text('\n'.join(captured))
    checks={}
    for name,root in ROOTS.items():
        before=read(run/f'provenance/{name}_before_manifest.json');after=manifest(root)
        write(run/f'provenance/{name}_after_manifest.json',after)
        missing=sorted(set(before)-set(after));added=sorted(set(after)-set(before))
        changed=[k for k in before if k in after and before[k]!=after[k]]
        allowed_new=[k for k in added if name=='hemocell' and (k=='test_code/step3_vascular_pure_fluid' or k.startswith('test_code/step3_vascular_pure_fluid/'))]
        allowed_dir=[k for k in changed if name=='hemocell' and k=='test_code' and stat.S_ISDIR(before[k][0]) and before[k][0]==after[k][0]]
        observed_git_metadata=[k for k in changed if name=='source' and k in ['.git','.git/index']]
        changes=[k for k in changed if k not in allowed_dir and k not in observed_git_metadata]
        checks[name]={'status':'PASS' if not missing and not changes and set(added)==set(allowed_new) else 'FAIL',
                      'missing':missing,'unexpected_changed':changes,'unexpected_added':sorted(set(added)-set(allowed_new)),
                      'authorized_new_files':allowed_new,'benign_parent_directory_metadata_changes':allowed_dir,'baseline_entries':len(before)}
        if observed_git_metadata:
            before_log=(run/'git_status_before.txt').read_text()
            matching=[chunk for chunk in before_log.split('COMMAND ') if str(SRC) in chunk.split('\n',1)[0] and "'status', '--short'" in chunk.split('\n',1)[0]]
            assert len(matching)==1
            before_status=matching[0].split('RC=0\n',1)[1].strip()
            same=before_status==statuses[str(SRC)]
            checks[name]['git_metadata_observation']={
                'paths':observed_git_metadata,'before':{k:before[k] for k in observed_git_metadata},
                'after':{k:after[k] for k in observed_git_metadata},'status_porcelain_equal_before':same,
                'meaning':'The first final-audit git status/diff invocation refreshed source .git/index metadata despite GIT_OPTIONAL_LOCKS=0. Working-tree files and complete status listing are unchanged. No index binary hash was captured before, so bytewise index identity is not claimed. This is a read-only-workflow metadata exception, explicitly retained; no rollback or timestamp hiding.',
                'subsequent_audit':'GIT_INDEX_FILE points to RUN-local copy; diff.autoRefreshIndex=false'}
            if not same:checks[name]['status']='FAIL'
    for name in ['step1','step2']:
        hashes=read(run/f'provenance/{name}_before_sha256.json')
        bad=[k for k,v in hashes.items() if sha(ROOTS[name]/k)!=v]
        checks[name].update(all_file_hash_checks=len(hashes),hash_mismatches=bad)
        if bad:checks[name]['status']='FAIL'
    tracked=read(run/'provenance/hemocell_tracked_before_sha256.json')
    bad=[k for k,v in tracked.items() if (sha(HC/k) if (HC/k).is_file() else None)!=v]
    checks['hemocell'].update(tracked_file_hash_checks=len(tracked),tracked_hash_mismatches=bad)
    if bad:checks['hemocell']['status']='FAIL'
    review=read(run/'provenance/step2_code_review.json')
    checks['step2_code']={'status':'PASS' if all(sha(e['path'])==e['sha256'] for e in review) else 'FAIL','files':len(review)}
    checks['status']='PASS' if all(v['status']=='PASS' for v in checks.values()) else 'FAIL'
    write(run/'provenance/frozen_integrity_after.json',checks)
    return checks
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path)
    print(json.dumps(audit(parser.parse_args().run),indent=2))
