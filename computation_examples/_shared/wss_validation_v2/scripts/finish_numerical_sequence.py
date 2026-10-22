"""Local ordered gates: finish two-mesh stage 3 only; pressure CFD SKIPPED_BY_USER.

Never substitutes a queued job, interpolation, or solver exit for CFD evidence.
Run alongside the single stage-three result watcher. Resume from recorded state.
"""
from pathlib import Path
import csv, hashlib, json, os, subprocess, sys, time, traceback

V = Path(__file__).resolve().parents[1]
R = '/workspace/wss_validation_v2_20260927T1230Z'
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
state = V/'logs/ordered_sequence_progress.json'

def save(**kw):
    state.write_text(json.dumps(dict(observed_unix=time.time(), **kw), indent=2)+'\n')

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run(args, log):
    with (V/'logs'/log).open('w') as f:
        subprocess.run(args, cwd=V, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)

def script(name, extra=(), log=None):
    run([sys.executable, '-B', str(V/'scripts'/name), *map(str, extra)], log or name+'.log')

def main():
    if not (V/'stage3/completion.json').exists():
        save(stage=3, status='WAITING_INDEPENDENT_BASELINE_AND_MEDIUM_RESULTS')
        while True:
            p = V/'logs/stage3_analysis_watch.json'
            d = read(p) if p.exists() else {}
            if 'FAILURE' in d.get('status', ''):
                raise RuntimeError('Stage-three watcher requires review: '+str(d))
            if d.get('status') == 'COMPLETE':
                break
            time.sleep(30)
        for name in ['vessel_baseline','vessel_medium']:
            case = V/'stage3'/name
            assert read(case/'reports/flow_quality.json')['accepted_final_and_log_checks']
            if name != 'vessel_baseline':
                ex = read(case/'reports/execution.json')
                assert ex['status']=='PASS' and ex['input_hashes_unchanged']
                assert read(case/'reports/linear_attempt_acceptance.json')['accepted']
        script('collect_results.py', log='ordered_collect_stage3.log')
        script('collect_geometry.py', log='ordered_geometry_stage3.log')
        script('reproduce_local_wall_evidence.py', log='ordered_local_subset_reproduction.log')
        comparisons = list(csv.DictReader((V/'data/vessel_mesh_comparison.csv').open()))
        targets = [r for r in comparisons if r['metric'] in ['mean_Pa','p05_Pa','p95_Pa']]
        assert len(targets)==33
        result = dict(status='ACTUAL_TWO_MESH_SENSITIVITY_ANALYSIS_COMPLETE',
                      completed_unix=time.time(), cases=['vessel_baseline','vessel_medium'],
                      reuse='Original H0 identity-checked; medium is actual new CFD; fine SKIPPED_BY_USER',
                      regional_metrics_compared=len(targets),
                      fine_status='SKIPPED_BY_USER',
                      reason='User cancelled fine vascular CFD due to time cost; three-grid criteria cannot be evaluated',
                      not_a_claim_of_mesh_independence=True,
                      comparison_csv_sha256=sha(V/'data/vessel_mesh_comparison.csv'))
        (V/'stage3/completion.json').write_text(json.dumps(result, indent=2)+'\n')
        script('make_figures.py', ['--part','vessel'], 'ordered_vessel_figures.log')

    # Stage four was explicitly cancelled by the user. This script contains no
    # case generation, SSH registration, solver launch, or pressure-result wait.
    cancellation=read(V/'stage4/user_cancellation.json')
    assert cancellation['status']=='SKIPPED_BY_USER'
    (V/'stage4/completion.json').write_text(json.dumps(cancellation,indent=2)+'\n')
    save(stage=3,status='AUTHORIZED_NUMERICAL_WORK_COMPLETE_FINAL_REPORT_REMAINS',fine='SKIPPED_BY_USER',pressure_sensitivity='SKIPPED_BY_USER')

if __name__=='__main__':
    try:main()
    except BaseException as e:
        save(status='REQUIRES_REVIEW',error=repr(e),traceback=traceback.format_exc())
        raise
