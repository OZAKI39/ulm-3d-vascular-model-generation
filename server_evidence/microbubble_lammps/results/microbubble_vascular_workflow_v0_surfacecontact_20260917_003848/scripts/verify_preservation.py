from pathlib import Path
import hashlib,json,subprocess,os,csv
S=Path(__file__).resolve().parents[1];B=S.parent/'bcflux_20260916_224835';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();records={}
D=S.parent/'wallslide_20260916_233923'
for rel,h in json.loads((S/'provenance/WALLSLIDE_IMMUTABLE_FILES.json').read_text()).items():assert sha(D/rel)==h;records[str(D/rel)]=h
for rel,h in json.loads((D/'provenance/BCFLUX_IMMUTABLE_FILES.json').read_text()).items():assert sha(B/rel)==h;records[str(B/rel)]=h
for name in ['DONOR_HASHES_BEFORE.json','ORIGINAL_HASHES_BEFORE.json']:
    for path,h in json.loads((B/'provenance'/name).read_text()).items():
        h=h['sha256'] if isinstance(h,dict) else h;assert sha(Path(path))==h;records[path]=h
for path,h in json.loads((S.parent/'adaptive_flux_20260916_222104/provenance/DONOR_HASHES_BEFORE.json').read_text()).items():
    h=h['sha256'] if isinstance(h,dict) else h;assert sha(Path(path))==h;records[path]=h
for e in json.loads((S/'validation/AUTHORITATIVE_INLET_FLOW_AUDIT.json').read_text())['chain']:assert sha(Path(e['source_file']))==e['sha256'];records[e['source_file']]=e['sha256']
before=json.loads((S/'provenance/GIT_BEFORE.json').read_text());after={}
for path in before:
    def g(args):return subprocess.check_output(['git','-C',path]+args,text=True,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'})
    after[path]={'branch':g(['branch','--show-current']),'commit':g(['rev-parse','HEAD']),'status':g(['status','--porcelain=v1','--untracked-files=all'])}
assert before==after
(S/'provenance/GIT_AFTER.json').write_text(json.dumps(after,indent=2)+'\n')
unchanged=['rigid_math.cpp','rigid_math.hpp','frozen_flow.cpp','frozen_flow.hpp','adaptive_injection.hpp','lifecycle.hpp','workflow_geometry.hpp']
for name in unchanged:assert sha(S/'src'/name)==sha(B/'src'/name)
for p in (D/'geometry').iterdir():
    if p.is_file():assert sha(S/'geometry'/p.name)==sha(p)
for directory in ['fields','inputs']:
    for p in (S/directory).iterdir():assert sha(p)==sha(B/directory/p.name)
(S/'validation/ORIGINAL_INPUT_PRESERVATION.json').write_text(json.dumps(dict(status='PASS',local_protected_files=len(records),bcflux_immutable_files=545,wallslide_immutable_files=292,production_geometry_hashes_unchanged=True,git_branch_commit_full_status_unchanged=True,unchanged_source_components=unchanged,source_modified_outside_new_stage=False,git_mutations=False),indent=2)+'\n')
with (S/'provenance/INPUT_MANIFEST.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['category','path','sha256','bytes'])
    for n,h in sorted(records.items()):w.writerow(['READ_ONLY_ORIGINAL_OR_DONOR',n,h,Path(n).stat().st_size])
    for dirname in ['geometry','fields','inputs','configs']:
        for p in sorted((S/dirname).iterdir()):
            if p.is_file():w.writerow(['STAGED_INPUT',str(p),sha(p),p.stat().st_size])
print('PRESERVATION_PASS',len(records))
