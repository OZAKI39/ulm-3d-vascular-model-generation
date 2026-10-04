"""Final audit of test identities, protected hashes and review deliverables."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/pfpd_at1_coupling_stability_audit_v1'


def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def cases(path):
    out={}
    for c in ET.parse(path).findall('.//testcase'):
        failure=c.find('failure'); error=c.find('error'); skipped=c.find('skipped')
        status='failed' if failure is not None else 'error' if error is not None else 'skipped' if skipped is not None else 'passed'
        out[c.attrib['classname'].replace('.','/')+'.py::'+c.attrib['name']]=dict(status=status,
            failure_message=None if failure is None else failure.attrib.get('message',''),
            failure_text=None if failure is None else failure.text)
    return out


def main():
    base=cases(OUT/'tests/BASELINE_CURRENT.xml');final=cases(OUT/'tests/FINAL_FULL.xml');diag=cases(OUT/'tests/DIAGNOSTIC.xml')
    changes=[k for k,v in base.items() if final.get(k)!=v]
    new=set(final)-set(base)
    summary=dict(baseline_passed=sum(v['status']=='passed' for v in base.values()),
        baseline_failed=sum(v['status']=='failed' for v in base.values()),
        final_passed=sum(v['status']=='passed' for v in final.values()),
        final_failed=sum(v['status']=='failed' for v in final.values()),
        final_skipped=sum(v['status']=='skipped' for v in final.values()),
        new_diagnostic_passed=sum(v['status']=='passed' for v in diag.values()),
        known_failed_IDs=[k for k,v in base.items() if v['status']=='failed'],
        exact_failure_messages_and_tracebacks_unchanged=not changes,
        changed_original_test_results=changes,
        all_new_tests_independently_passed=new==set(diag) and all(v['status']=='passed' for v in diag.values()),
        new_regressions=[k for k in new if final[k]['status']!='passed'])
    write(OUT/'tests/COMPARISON.json',summary)
    manifest=read(OUT/'provenance/PRE_RUN_MANIFEST.json');repo=Path(manifest['repo_root'])
    modified=[];missing=[];post={};link_changes=[]
    for name,r in manifest['files'].items():
        p=repo/name
        if not p.is_file() or p.is_symlink():missing.append(name);continue
        value=dict(sha256=sha(p),bytes=p.stat().st_size);post[name]=value
        if value!=r:modified.append(name)
    for name,target in manifest['symlinks'].items():
        p=repo/name
        if not p.is_symlink() or os.readlink(p)!=target:link_changes.append(name)
    new_managed={}
    excluded=[]
    own_final={OUT/'provenance/POST_RUN_MANIFEST.json',OUT/'provenance/PRESERVATION_CHECK.json',OUT/'provenance/FINAL_REVIEW.json'}
    for folder,dirs,files in os.walk(repo):
        dirs[:]=[d for d in dirs if d!='.git' and not (Path(folder)/d).is_symlink()]
        here=Path(folder)
        if here==OUT/'provenance/runtime_cache' or (here.parent==OUT/'tests' and here.name.endswith('_tmp')):
            excluded.append(str(here.relative_to(repo)));dirs[:]=[];continue
        for name in files:
            p=here/name;rel=str(p.relative_to(repo))
            if rel in manifest['files'] or p.is_symlink() or p in own_final:continue
            new_managed[rel]=dict(sha256=sha(p),bytes=p.stat().st_size)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    preservation=dict(status='PASSED' if not(modified or missing or link_changes) and head==manifest['head'] else 'FAILED',
        original_regular_files_checked=len(manifest['files']),original_symlinks_checked=len(manifest['symlinks']),
        existing_files_modified=len(modified),existing_files_missing=len(missing),symlinks_changed=len(link_changes),
        production_files_modified=sum(n.startswith('clot_clearance/pd_clot/') for n in modified),
        historical_result_files_modified=sum(n.startswith(('clot_clearance/results/','clot_clearance/verification/')) for n in modified),
        modified=modified,missing=missing,changed_symlinks=link_changes,head_before=manifest['head'],head_after=head)
    write(OUT/'provenance/PRESERVATION_CHECK.json',preservation)
    write(OUT/'provenance/POST_RUN_MANIFEST.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        repo_root=str(repo),head=head,original_files=post,original_symlinks=manifest['symlinks'],new_managed_files=new_managed,
        excluded_new_runtime_directories=excluded,self_referential_files_excluded=[str(p.relative_to(repo)) for p in sorted(own_final)]))
    names=['iteration_amplification','fixed_point_spectrum','observed_vs_eigenmode','monolithic_vs_staggered_energy',
           'smallest_hessian_eigenvalue','unstable_mode_map','substep_diagnostic','tension_gate_ordering']
    deliverables=all((OUT/'figures'/f'{name}.{ext}').stat().st_size>1000 for name in names for ext in ['png','pdf'])
    report=(OUT/'REVIEW_COUPLING_STABILITY_AUDIT.md').is_file()
    passed=preservation['status']=='PASSED' and not changes and summary['all_new_tests_independently_passed'] and not summary['new_regressions'] and deliverables and report
    write(OUT/'provenance/FINAL_REVIEW.json',dict(passed=passed,preservation=preservation['status'],tests=summary,
        figures_png_pdf_pairs=8 if deliverables else None,report_present=report,
        status=read(OUT/'STATUS.json')))
    print(json.dumps(dict(passed=passed,preservation=preservation,tests=summary,figures_png_pdf_pairs=8 if deliverables else None),indent=2))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
