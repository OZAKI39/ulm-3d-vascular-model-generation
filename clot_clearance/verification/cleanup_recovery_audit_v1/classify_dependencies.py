"""Classify the frozen 123-item population against unchanged test contracts."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

A = Path(__file__).resolve().parent
PROJECT = A.parents[1]
ROOT = A.parents[2]
OLD = A.parent / 'repository_cleanup_v1'
ROWS = json.loads((OLD / 'RECOVERY_UNAVAILABLE.json').read_text())
SEMANTIC = {r['path']: r for r in json.loads((A / 'SEMANTIC_RECOVERY_COMPARISON.json').read_text())}
RULES = [
    ('pfpd_at1_spatial_discretization_audit_v2/provenance/PRE_RUN_MANIFEST.json',
     'tests/test_spatial_discretization_v2.py::test_legacy_files_are_unchanged', 63,
     ('results/pfpd_at1_evolution_v1/', 'verification/pfpd_at1_branch_stability_tracking_v1/', 'verification/pfpd_at1_spatial_discretization_audit_v1/')),
    ('pfpd_at1_spatial_discretization_audit_v2/provenance/PRE_RUN_MANIFEST.json',
     'tests/test_spatial_discretization_v2.py::test_historical_sources_and_tests_unchanged', 114,
     ('pd_clot/', 'tests/', 'scripts/')),
    ('pfpd_at1_spectral_split_v1/provenance/PRE_RUN_MANIFEST.json',
     'tests/test_pfpd_spectral_domain.py::test_all_preexisting_production_and_candidate_sources_unchanged', 29,
     ('pd_clot/', 'scripts/', 'tests/')),
    ('pfpd_at1_c5_stability_v1/PRE_RUN_MANIFEST.json',
     'tests/test_c5_stability_bracket.py::test_saved_roots_production_and_existing_tests_unchanged', 34,
     ('pd_clot/', 'tests/', 'verification/pfpd_at1_petsc_vi_solver_v1/continuation/')),
    ('pfpd_at1_c7_early_branch_v1/PRE_RUN_MANIFEST.json',
     'tests/test_c7_early_branch.py::test_production_and_historical_geometry_unchanged', 40,
     ('pd_clot/',)),
]

def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args], text=True,
                                   env=dict(os.environ, GIT_OPTIONAL_LOCKS='0', GIT_NO_LAZY_FETCH='1'))

def main():
    tracked = set(git('ls-files').splitlines())
    history = git('log', '--all', '--reflog', '--name-only', '--format=')
    (A / 'GIT_EVER_TRACKED_PATHS.log').write_text(history)
    ever = set(history.splitlines())
    # Search exact path and recognizable fixture-relative path in first-party text.
    # Manifest membership is handled separately; generic cache access is recorded below.
    needles = {}
    for row in ROWS:
        rel = row['path'].removeprefix('clot_clearance/')
        variants = [rel, row['path']]
        if '/tests/' in rel:
            parts = rel.split('/tests/', 1)[1].split('/', 1)
            if len(parts) > 1:
                variants.append(parts[1])
        needles[row['path']] = variants
    references = {r['path']: dict(test=[], code=[], review=[]) for r in ROWS}
    scanned = []
    for folder, dirs, files in os.walk(PROJECT, followlinks=False):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__', 'external', 'source_snapshot',
                    'runtime_cache', 'baseline_tmp', 'final_tmp', 'repository_cleanup_v1',
                    'cleanup_recovery_audit_v1') and not d.startswith('build_')]
        for name in files:
            p = Path(folder) / name
            rel = str(p.relative_to(PROJECT))
            if p.is_symlink() or 'MANIFEST' in name.upper():
                continue
            role = ('test' if rel.startswith('tests/') and p.suffix == '.py' else
                    'code' if rel.startswith(('pd_clot/', 'scripts/', 'vendor/')) and p.suffix == '.py' else
                    'review' if p.suffix.lower() == '.md' or any(t in name.upper() for t in ('SUMMARY', 'STATUS', 'REVIEW')) else None)
            if role is None or p.stat().st_size > 6_000_000:
                continue
            text = p.read_text(errors='replace')
            scanned.append(rel)
            for row in ROWS:
                for needle in needles[row['path']]:
                    if needle in text:
                        references[row['path']][role].append(rel)
                        break
    manifests = {name: json.loads((PROJECT / 'verification' / name).read_text())['files'] for name, *_ in RULES}
    out, evidence = [], []
    for row in ROWS:
        rel = row['path'].removeprefix('clot_clearance/')
        deps = []
        for manifest, test, line, prefixes in RULES:
            if rel.startswith(prefixes) and rel in manifests[manifest]:
                assert manifests[manifest][rel] == row['hash'], rel
                deps.append(dict(test=test, source_line=line, manifest='verification/' + manifest,
                                 operation='read_bytes() then SHA256 equality', expected_sha256=row['hash']))
        assert deps, ('unclassified test dependency', rel)
        exact = (ROOT / row['path']).is_file() and hashlib.sha256((ROOT / row['path']).read_bytes()).hexdigest() == row['hash']
        semantic = SEMANTIC.get(row['path'], {})
        suffix = Path(row['path']).suffix
        secondary = 'REGENERABLE_CACHE' if suffix in ('.nbc', '.nbi') or 'fontlist' in rel else 'TEMPORARY_ARTIFACT'
        code = references[row['path']]['code']
        if suffix in ('.nbc', '.nbi'):
            code = code + ['OPTIONAL_NUMBA_CACHE: @njit(cache=True); runtime regenerates cache, original bytes are not scientifically required']
        scientific_role = ('optional compiler/font cache; no accepted-state data' if secondary == 'REGENERABLE_CACHE' else
                           'historical pytest fixture artifact; scientific run records remain at same fixture; not accepted C5/C7 state')
        status = 'EXACT_BYTE_RECOVERED' if exact else 'CRITICAL_UNRECOVERED'
        reconstruction = ('EXACT_VERIFIED' if exact else
                          'SCIENTIFIC_EQUIVALENT_STAGED; original hash test unresolved' if semantic.get('classification') == 'SCIENTIFIC_EQUIVALENT' else
                          'PARTIAL_NOT_EQUIVALENT' if semantic.get('classification') == 'NOT_EQUIVALENT' else
                          'CACHE_REGENERABLE_BUT_ORIGINAL_SERIALIZATION_NOT_VERIFIED')
        out.append(dict(path=row['path'], expected_sha256=row['hash'], expected_size=row['size_bytes'],
            source_manifest=';'.join(sorted({d['manifest'] for d in deps})), git_tracked=row['path'] in tracked,
            extension=suffix, referenced_by_test=';'.join(d['test'] for d in deps),
            referenced_by_review=';'.join(references[row['path']]['review']) or 'NO_LITERAL_REFERENCE_FOUND',
            referenced_by_code=';'.join(code) or 'NO_LITERAL_RUNTIME_REFERENCE_FOUND',
            scientific_role=scientific_role, reconstructable=reconstruction, recovery_status=status,
            dependency_class='DIRECT_TEST_DEPENDENCY', secondary_class=secondary,
            git_ever_tracked=row['path'] in ever,
            semantic_status=semantic.get('classification', 'NOT_NEEDED_EXACT_RESTORED' if exact else 'NOT_RECONSTRUCTED')))
        evidence.append(dict(path=row['path'], direct_test_dependencies=deps,
                             literal_references=references[row['path']], current_scientific_dependency='none found',
                             classification_reason='The test itself reads artifact bytes through a manifest. Manifest indirection does not make this MANIFEST_ONLY.',
                             semantic_recovery=semantic))
    with (A / 'RECOVERY_MISSING_123.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    critical = [r for r in out if r['recovery_status'] == 'CRITICAL_UNRECOVERED']
    with (A / 'CRITICAL_UNRECOVERED.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(critical)
    (A / 'DEPENDENCY_EVIDENCE.json').write_text(json.dumps(dict(
        scan_files=scanned, scope='First-party tests/runtime/review/status/summary; exact and fixture-relative literal paths plus explicit manifest-loop resolution.',
        limitations='Dynamic reads outside inspected code cannot be disproven by string search; optional compiler cache loading is recorded explicitly.',
        items=evidence), indent=2) + '\n')
    exceptions = [dict(original_path=r['path'], expected_hash=r['expected_sha256'],
                       reason_lost='Removed by previous repository_cleanup_v1 after indirect test-manifest dependency was missed.',
                       scientific_dependency='none for current accepted C5/C7 results or production runtime',
                       historical_scientific_information_loss=('Some historical fixture extrema are not recoverable' if r['semantic_status']=='NOT_EQUIVALENT' else 'see semantic recovery evidence'),
                       actual_test_dependency=r['referenced_by_test'], exception_authorizes_test_change=False,
                       current_replacement=SEMANTIC.get(r['path'], {}).get('staged'),
                       semantic_recovery_status=r['semantic_status'], disposition='CRITICAL_UNRECOVERED: old test byte contract unresolved') for r in critical]
    (A / 'RECOVERY_EXCEPTION_REPORT.json').write_text(json.dumps(dict(
        status='NOT_A_WAIVER', old_manifests_unchanged=True, tests_unchanged=True,
        note='These artifacts are scientifically unnecessary for current accepted states, but remain actual old-test dependencies. No no-recovery-needed exception is approved.',
        items=exceptions), indent=2) + '\n')
    from collections import Counter
    summary=dict(population=len(out), classified=len(out), unknown=0,
                 statuses=dict(Counter(r['recovery_status'] for r in out)),
                 primary_classes=dict(Counter(r['dependency_class'] for r in out)),
                 secondary_classes=dict(Counter(r['secondary_class'] for r in out)),
                 current_git_tracked=sum(r['git_tracked'] for r in out), ever_git_tracked=sum(r['git_ever_tracked'] for r in out),
                 critical_reason='Direct old-test SHA contracts still unsatisfied; scientific-equivalent staged JSON is not an exact or full-semantic replacement for that contract.')
    (A / 'CLASSIFICATION_SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary), flush=True)

if __name__ == '__main__':
    main()
