"""Review original test identities, historical hashes and new branch artifacts."""
from pathlib import Path
import os
import subprocess
import json
from datetime import datetime,timezone
import numpy as np
from scripts.review_coupling_stability_audit import cases,sha,read,write

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/pfpd_at1_branch_stability_tracking_v1'


def main():
    baseline=cases(OUT/'tests/BASELINE_CURRENT.xml');final=cases(OUT/'tests/FINAL_FULL.xml');new=cases(OUT/'tests/DIAGNOSTIC_FINAL.xml')
    changes=[k for k,v in baseline.items() if final.get(k)!=v]
    test_summary=dict(baseline_passed=sum(v['status']=='passed' for v in baseline.values()),
        known_failed_IDs=[k for k,v in baseline.items() if v['status']=='failed'],
        final_passed=sum(v['status']=='passed' for v in final.values()),final_failed=sum(v['status']=='failed' for v in final.values()),
        final_skipped=sum(v['status']=='skipped' for v in final.values()),new_diagnostic_passed=sum(v['status']=='passed' for v in new.values()),
        original_test_changes=changes,exact_failure_messages_and_tracebacks_unchanged=not changes,
        all_new_tests_pass=set(final)-set(baseline)==set(new) and all(v['status']=='passed' for v in new.values()),
        preliminary_new_test_failure_retained='DIAGNOSTIC.log / DIAGNOSTIC.xml; finite-difference cancellation, no original test or numerical solver tolerance changed')
    write(OUT/'tests/COMPARISON.json',test_summary)
    manifest=read(OUT/'provenance/PRE_RUN_MANIFEST.json');repo=Path(manifest['repo_root']);modified=[];missing=[];links=[];post={}
    for name,value in manifest['files'].items():
        p=repo/name
        if not p.is_file() or p.is_symlink():missing.append(name);continue
        record=dict(sha256=sha(p),bytes=p.stat().st_size);post[name]=record
        if record!=value:modified.append(name)
    for name,target in manifest['symlinks'].items():
        p=repo/name
        if not p.is_symlink() or os.readlink(p)!=target:links.append(name)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    preservation=dict(status='PASSED' if not(modified or missing or links) and head==manifest['head'] else 'FAILED',
        original_files_checked=len(manifest['files']),original_symlinks_checked=len(manifest['symlinks']),
        existing_files_modified=len(modified),production_files_modified=sum(n.startswith('clot_clearance/pd_clot/') for n in modified),
        historical_files_modified=sum(n.startswith(('clot_clearance/results/','clot_clearance/verification/')) for n in modified),
        existing_tests_modified=sum(n.startswith('clot_clearance/tests/') for n in modified),
        modified=modified,missing=missing,symlinks_changed=links,head=head)
    write(OUT/'provenance/PRESERVATION_CHECK.json',preservation)
    closure={}
    for k in ['A','B']:
        b=read(OUT/f'critical_bracket/SCHEDULE_{k}.json')
        with np.load(OUT/'branch_states'/b['lower']['state']/'STATE.npz') as low,np.load(OUT/'branch_states'/b['upper']['state']/'STATE.npz') as high:
            closure[k]=bool(np.array_equal(low['s'],high['lower']))
    write(OUT/'comparisons/FINAL_BRACKET_BOUND_COHERENCE.json',dict(upper_physical_lower_equals_final_lower_phase=closure))
    names=['stability_eigenvalue_vs_strain','critical_bracket_zoom','branch_force_phase','active_set_evolution',
        'curvature_along_critical_mode','critical_mode_map','stationarity_quality','schedule_consistency']
    figures=all((OUT/'figures'/f'{n}.{ext}').is_file() and (OUT/'figures'/f'{n}.{ext}').stat().st_size>1000 for n in names for ext in ['png','pdf'])
    report=(OUT/'REVIEW_BRANCH_STABILITY_TRACKING.md').is_file()
    declared=sha(OUT/'provenance/PREDECLARED_CONFIG.json')==(OUT/'provenance/PREDECLARED_CONFIG.sha256').read_text().strip()
    passed=preservation['status']=='PASSED' and not changes and test_summary['all_new_tests_pass'] and all(closure.values()) and figures and report and declared
    write(OUT/'provenance/FINAL_REVIEW.json',dict(passed=passed,tests=test_summary,preservation=preservation['status'],
        final_bracket_bound_coherence=closure,predeclared_config_unchanged=declared,figure_pairs=8 if figures else None,report_present=report,status=read(OUT/'STATUS.json')))
    managed={};excluded=[];self_files={OUT/'provenance/POST_RUN_MANIFEST.json'}
    for folder,dirs,names in os.walk(repo):
        here=Path(folder);dirs[:]=[d for d in dirs if d!='.git' and not (here/d).is_symlink()]
        if here==OUT/'provenance/runtime_cache' or (here.parent==OUT/'tests' and here.name.endswith('_tmp')):
            dirs[:]=[];excluded.append(str(here.relative_to(repo)));continue
        for n in names:
            p=here/n;relative=str(p.relative_to(repo))
            if relative in manifest['files'] or p.is_symlink() or p in self_files:continue
            managed[relative]=dict(sha256=sha(p),bytes=p.stat().st_size)
    write(OUT/'provenance/POST_RUN_MANIFEST.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),repo_root=str(repo),head=head,
        original_files=post,original_symlinks=manifest['symlinks'],new_managed_files=managed,excluded_new_runtime_directories=excluded,
        self_referential_file_excluded=str((OUT/'provenance/POST_RUN_MANIFEST.json').relative_to(repo))))
    print(json.dumps(dict(passed=passed,preservation=preservation,tests=test_summary,bracket_bound_coherence=closure,figure_pairs=8 if figures else None),indent=2))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
