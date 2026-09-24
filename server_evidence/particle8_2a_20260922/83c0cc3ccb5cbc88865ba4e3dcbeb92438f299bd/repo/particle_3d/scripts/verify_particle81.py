#!/usr/bin/env python3
"""Permanent source/evidence-bound P8.1 tests and relevant upstream regressions."""
from pathlib import Path
import argparse,subprocess,sys,xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle81_simulation import OUTPUT,dump,lock_upstream
from particle_3d.particle8_replay import REPO,digest,read
UPSTREAM=['particle1','particle6_5','particle7','particle8','particle8_full3d']


def source_files():
    paths=[]
    for pattern in ['PARTICLE8_1_README.md','scripts/*particle81*.py','src/particle_3d/particle81*.py','tests/particle81/*.py']:
        paths.extend((REPO/'particle_3d').glob(pattern))
    return {str(p.relative_to(REPO)):digest(p) for p in sorted(paths)}


def artifact_files():
    # Diagnostic archives are inventoried separately; they are not the formal cohort.
    names=[]
    for folder in ['data','trajectories','frames','animations','keyframes','figures','inspection']:
        names.extend(p for p in (OUTPUT/folder).glob('*') if p.is_file() and p.name not in
                     ['test_run.json','visual_inspection.json','render_determinism.json'])
    for relative in ['data/reference_original_run/trajectories','data/cache_v2_benchmark/trajectories','data/sensitivity']:
        names.extend(p for p in (OUTPUT/relative).rglob('*') if p.is_file())
    return {str(p.relative_to(OUTPUT)):digest(p) for p in sorted(names)}


def run_suite(suite):
    xml=OUTPUT/'logs'/f'tests_{suite}.xml';log=OUTPUT/'logs'/f'tests_{suite}.log'
    cmd=[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',f'--junitxml={xml}',f'particle_3d/tests/{suite}']
    with log.open('w') as stream:result=subprocess.run(cmd,cwd=REPO,stdout=stream,stderr=subprocess.STDOUT)
    suites=ET.parse(xml).getroot().findall('testsuite')
    counts={key:sum(int(s.get(key,0)) for s in suites) for key in ['tests','failures','errors','skipped']}
    record=dict(suite=suite,command=cmd,exit_code=result.returncode,**counts,
        junit=str(xml.relative_to(OUTPUT)),log=str(log.relative_to(OUTPUT)),junit_sha256=digest(xml),log_sha256=digest(log))
    print(suite,counts,flush=True);return record


def passed(results):
    return all(r['exit_code']==0 and r['tests']>0 and not (r['failures']+r['errors']+r['skipped']) for r in results)


def upstream_results(refresh=False):
    lock_upstream();lock=read(OUTPUT/'data/upstream_lock.json');receipt=OUTPUT/'data/upstream_test_run.json'
    if receipt.exists() and not refresh:
        data=read(receipt)
        valid=data['scope_sha256']==lock['sha256'] and passed(data['suites']) and [r['suite'] for r in data['suites']]==UPSTREAM
        valid &= all(digest(OUTPUT/r[key])==r[key+'_sha256'] for r in data['suites'] for key in ['junit','log'])
        if valid:return data['suites']
    results=[run_suite(s) for s in UPSTREAM];lock_upstream()
    dump(receipt,dict(suites=results,all_pass=passed(results),scope_sha256=lock['sha256'],
        reuse_rule='Run during this P8.1 development; reusable only while every locked upstream source/test/input and test receipt is byte-identical'))
    return results


def main(upstream_only=False):
    if upstream_only:return 0 if passed(upstream_results(refresh=True)) else 1
    # Upstream regressions can run while independent MB workers are integrating.
    # Their actual current-stage executions are reused only under the full lock.
    upstream=upstream_results();before=source_files();artifacts=artifact_files()
    results=[run_suite('particle81'),*upstream]
    unchanged=before==source_files() and artifacts==artifact_files()
    ok=unchanged and passed(results)
    dump(OUTPUT/'data/test_run.json',dict(suites=results,all_pass=ok,inputs_unchanged_during_tests=unchanged,
        source_sha256=before,artifact_sha256=artifacts,upstream_preserved_files=lock_upstream(),
        upstream_execution_receipt='data/upstream_test_run.json'))
    return 0 if ok else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--upstream-only',action='store_true');args=p.parse_args()
    sys.exit(main(args.upstream_only))
