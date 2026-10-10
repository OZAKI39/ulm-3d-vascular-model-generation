"""Finalize evidence only; never launch mechanics or modify old files."""
from pathlib import Path
import sys,json,hashlib,shutil,zipfile,subprocess,time
import psutil
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from pd_clot.lammps_fragmentation.targeted_contact_transfer import O,R,PREVIOUS,OLD,GAUSSIAN100
from pd_clot.lammps_fragmentation.config import ROOT,dump

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda:f.read(8*1024*1024),b''):h.update(data)
    return h.hexdigest()
start=time.perf_counter();assert not (O/'evidence.zip').exists()
for filename in ['checks/ROUND_AUDIT.json','checks/HISTORY_INTEGRITY.json','VISUALIZATION_QA.json']:
    assert json.loads((O/filename).read_text())['status']=='PASS'
status=json.loads((O/'FINAL_STATUS.json').read_text());assert status['NEW_FORMAL_SOLVES']==4
new={'detached_fragment_metrics.py','targeted_contact_transfer.py','targeted_map_refresh.py','contact_transfer_visualization.py'};code=[]
for path in sorted((ROOT/'pd_clot/lammps_fragmentation').glob('*.py')):
    dest=O/'source_archive/python/pd_clot/lammps_fragmentation'/path.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,dest)
    code.append(dict(original_path=str(path),archived_path=str(dest.relative_to(O)),sha256=sha(path),role='NEW_THIS_ROUND' if path.name in new else 'REUSED_UNMODIFIED'))
for directory,filenames in [(OLD,['INTERFACE_GEOMETRY.npz','INTERFACE_GEOMETRY.json','CONFIGURATION_FROZEN.json','SOURCE_REGISTER.json','COHERENT_FAILURE_CARD_LIU2024_V1.json','ELASTIC_MAPPING_STATUS.json','REVIEW_LAMMPS_COHERENT_FAILURE_FORCE_TIME_MAP_V1.md']),(PREVIOUS,['BUILD_VARIANT.json','LAMMPS_CONTACT_PATCH.diff','OFFICIAL_SOURCE_CHECK.md','checks/CONTACT_STYLE_SMOKE.json','REVIEW_LAMMPS_CONTACT_CASCADE_DIAGNOSIS_AND_MINIMAL_FIX_V1.md'])]:
    for filename in filenames:
        dest=O/'source_archive'/directory.name/filename;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(directory/filename,dest)
dump(O/'SOURCE_ARCHIVE.json',dict(code=code,new_modules=sorted(new),new_packages_installed=[],new_LAMMPS_build=False))
historical=json.loads((O/'HISTORY_PROTECTION.json').read_text());groups={}
for row in historical:
    p=Path(row['path']);rel=p.relative_to(ROOT.parent);group=str(Path(*rel.parts[:3])) if rel.parts[0]=='clot_clearance' and rel.parts[1] in ['verification','results'] else 'other_protected_sources_and_executables'
    groups.setdefault(group,[]).append(row)
digests={}
for group,rows in groups.items():
    canonical=''.join(f"{r['path']}\t{r['bytes']}\t{r['sha256']}\n" for r in sorted(rows,key=lambda r:r['path']))
    digests[group]=dict(files=len(rows),bytes=sum(r['bytes'] for r in rows),starting_file_manifest_digest_sha256=hashlib.sha256(canonical.encode()).hexdigest(),end_file_by_file_hash_match=True)
dump(O/'HISTORICAL_DIRECTORY_HASHES.json',dict(groups=digests,definition='Aggregate of sorted path, bytes and SHA256 in the starting protection list; ending file hashes independently checked by audit_results.py'))
(O/'ENDING_GIT_STATE.md').write_text('# Ending workspace state\n\nNo commit, push, PR, cleanup, reset or stash. Existing changes retained.\n\n```text\n'+subprocess.check_output(['git','branch','--show-current'],cwd=ROOT.parent,text=True)+subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT.parent,text=True)+subprocess.check_output(['git','status','--short'],cwd=ROOT.parent,text=True)+'```\n')
active=[]
for proc in psutil.process_iter(['name','pid']):
    if proc.info['name'] in ['lmp','mpirun','mpiexec']:
        try:
            if str(O) in proc.cwd():active.append(proc.info)
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
assert not active
dump(O/'checks/STOP_STATE.json',dict(status='STOPPED_FOR_HUMAN_REVIEW',active_round_solvers=active,new_formal_cases=4,remaining_authorized_automatic_solves=0))
exclude={'evidence.zip','evidence.tmp.zip','evidence.zip.sha256','HASH_MANIFEST.json','DELIVERY_CHECK.json','package.log'}
large={'trajectory.lammpstrj','surviving_bonds.local','final.restart','components.csv'}
paths=sorted({p for root in [O,R] for p in root.rglob('*') if p.is_file() and p.name not in exclude and '__pycache__' not in p.parts})
# One historical saved-state movie reference; other old native raw files stay referenced.
paths+=sorted(p for p in GAUSSIAN100.iterdir() if p.is_file())
summaries=json.loads((O/'ALL_CASE_SUMMARIES.json').read_text());referenced=[]
for q in summaries.values():
    p=Path(q['source_path'])
    if p==GAUSSIAN100 or O in p.parents:continue
    referenced.extend(sorted(x for x in p.iterdir() if x.is_file()))
rows=[];members=[]
for path in sorted(set(paths+referenced)):
    included=path not in referenced and path.name not in large
    rel=str(path.relative_to(ROOT));row=dict(path=str(path),relative_path=rel,bytes=path.stat().st_size,sha256=sha(path),in_evidence_zip=included,retained_locally=True)
    rows.append(row)
    if included:members.append((path,rel))
required=['TASK_REQUEST.md','STARTING_GIT_STATE.md','FROZEN_CONFIGURATION.json','SOURCE_REGISTER.md','SIZE_THRESHOLD_POLICY.md','CONTACT_VARIANT_PROVENANCE.md','HISTORICAL_GAUSSIAN_REANALYSIS.json','PER_FRAME_DETACHED_COMPONENTS.csv','FRAGMENT_LIFETIMES.csv','FRAGMENT_SIZE_SUMMARY.csv','CROSS_PATTERN_CONTACT_DECISION.json','CONDITIONAL_CASE_DECISION.json','TARGETED_OLD_NEW_MAP_COMPARISON.csv','PERFORMANCE_REAL_CASES.json','FINAL_STATUS.json','REVIEW_LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1.md']
assert all((O/f).is_file() for f in required)
dump(O/'HASH_MANIFEST.json',dict(hash='SHA256',files=rows,large_raw_and_redundant_component_tables_retained=True,previous_complete_tree_hashes='HISTORY_PROTECTION.json',manifest_and_delivery_self_excluded=True))
members.append((O/'HASH_MANIFEST.json',str((O/'HASH_MANIFEST.json').relative_to(ROOT))))
with zipfile.ZipFile(O/'evidence.tmp.zip','w',zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
    for p,rel in members:z.write(p,rel)
(O/'evidence.tmp.zip').rename(O/'evidence.zip');checksum=sha(O/'evidence.zip');(O/'evidence.zip.sha256').write_text(checksum+'  evidence.zip\n')
with zipfile.ZipFile(O/'evidence.zip') as z:
    assert z.testzip() is None
    for row in rows:
        if row['in_evidence_zip']:assert hashlib.sha256(z.read(row['relative_path'])).hexdigest()==row['sha256']
    count=len(z.namelist())
for row in rows:assert sha(Path(row['path']))==row['sha256'],row['path']
dump(O/'DELIVERY_CHECK.json',dict(status='PASS',required_files_present=True,zip_CRC='PASS',member_SHA256='PASS',disk_manifest_recheck='PASS',history_1312_files='PASS',media_QA='PASS',scientific_audit='PASS',formal_solves=4,zip_members=count,zip_bytes=(O/'evidence.zip').stat().st_size,zip_sha256=checksum,raw_omitted_and_retained=sum(not r['in_evidence_zip'] for r in rows),active_round_solvers=[],stop='COMPLETE_STOPPED_FOR_HUMAN_REVIEW',packaging_and_validation_wall_s=time.perf_counter()-start))
print('DELIVERY_PASS',count,'members',round((O/'evidence.zip').stat().st_size/1024**2,2),'MiB',checksum,flush=True)
