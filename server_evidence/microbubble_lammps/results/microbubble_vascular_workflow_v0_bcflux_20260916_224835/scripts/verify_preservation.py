from pathlib import Path
import hashlib,json,subprocess,os,csv
S=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
groups={};all_records={}
for name in ['DONOR_HASHES_BEFORE.json','ORIGINAL_HASHES_BEFORE.json']:
    d=json.loads((S/'provenance'/name).read_text());groups[name]=len(d)
    for n,h in d.items():
        expected=h['sha256'] if isinstance(h,dict) else h;p=Path(n);assert sha(p)==expected,n;all_records[n]=expected
old=S.parent/'adaptive_flux_20260916_222104'
d=json.loads((old/'provenance/DONOR_HASHES_BEFORE.json').read_text());groups['preceding_planar_stage']=len(d)
for n,h in d.items():
    expected=h['sha256'] if isinstance(h,dict) else h
    assert sha(Path(n))==expected,n;all_records[n]=expected
audit=json.loads((S/'validation/AUTHORITATIVE_INLET_FLOW_AUDIT.json').read_text());groups['authoritative_BC_source_files']=len(audit['chain'])
for e in audit['chain']:assert sha(Path(e['source_file']))==e['sha256'];all_records[e['source_file']]=e['sha256']
gitbefore=json.loads((S/'provenance/GIT_BEFORE.json').read_text());after={}
for repo in gitbefore:
    def git(*a):return subprocess.check_output(['git','-C',repo,*a],env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'},text=True)
    after[repo]=dict(branch=git('branch','--show-current'),commit=git('rev-parse','HEAD'),status=git('status','--porcelain=v1','--untracked-files=all'))
assert after==gitbefore
(S/'provenance/GIT_AFTER.json').write_text(json.dumps(after,indent=2)+'\n')
(S/'validation/ORIGINAL_INPUT_PRESERVATION.json').write_text(json.dumps(dict(status='PASS',groups=groups,unique_read_only_files=len(all_records),git_branch_commit_full_status_unchanged=True,git_workspaces=list(after),no_git_mutating_commands=True),indent=2)+'\n')
with (S/'provenance/INPUT_MANIFEST.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['category','path','sha256','bytes'])
    for n,h in sorted(all_records.items()):w.writerow(['READ_ONLY_ORIGINAL_OR_DONOR',n,h,Path(n).stat().st_size])
    for dirname in ['geometry','fields','inputs']:
        for p in sorted((S/dirname).iterdir()):
            if p.is_file():w.writerow(['STAGED_INPUT',str(p),sha(p),p.stat().st_size])
for n in ['rigid_math.cpp','rigid_math.hpp','frozen_flow.cpp','frozen_flow.hpp','wall_distance.cpp','wall_distance.hpp']:
    assert sha(S/'src'/n)==sha(S/'provenance/stable_rigid_source'/n)
print('ORIGINAL_INPUT_AND_GIT_PRESERVATION_PASS',groups)
