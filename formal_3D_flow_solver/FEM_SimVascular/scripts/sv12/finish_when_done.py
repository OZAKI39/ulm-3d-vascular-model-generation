#!/usr/bin/env python3
"""Observe live nonlinear history, then finish all authorized SV1.2 evidence.

This process never launches, stops, or reconfigures the native solver. The
production runner owns the frozen run, safety stops and optional extension.
"""
import csv
import json
import os
import subprocess
import sys
import time
import traceback
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from sv_validation.sv12 import REPORT, OUTPUT, LOG, load
from sv_validation.provenance import write_json, now

def main():
    positions, pending, nonlinear = {}, {}, {}
    last_print = 0
    while not (REPORT/'flow_execution.json').exists():
        changed = False
        for name in ('production', 'extension'):
            path = OUTPUT/'qc'/(name+'_linear_history.jsonl')
            if not path.exists():
                continue
            with path.open() as stream:
                stream.seek(positions.get(name, 0))
                while True:
                    offset = stream.tell()
                    line = stream.readline()
                    if not line or not line.endswith('\n'):
                        positions[name] = offset
                        break
                    row = json.loads(line)
                    step = row['step']
                    initial = row['petsc_monitor'][0]['residual_norm'] if row.get('petsc_monitor') else None
                    if step not in pending:
                        pending[step] = initial
                    ratios = (row['nonlinear_Ri_over_R0'], row['nonlinear_Ri_over_R1'])
                    converged = row['nonlinear_iteration'] >= 2 and all(v is not None for v in ratios) and min(ratios) <= 1e-10
                    nonlinear[step] = {
                        'step': step, 'time_s': row['time_s'], 'iterations': row['nonlinear_iteration'],
                        'initial_residual': pending[step], 'final_residual': initial,
                        'Ri_over_R0': ratios[0], 'Ri_over_R1': ratios[1],
                        'converged': converged, 'source': name,
                        'observation': 'Residual before current Newton correction, as reported in native/PETSc history'}
                    changed = True
        if changed:
            rows = list(nonlinear.values())
            with (OUTPUT/'qc/nonlinear_live.csv').open('w') as out:
                writer = csv.DictWriter(out, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            write_json(REPORT/'nonlinear_live.json', {
                'timestamp': now(), 'observed_steps': len(rows),
                'completed_converged_steps': sum(r['converged'] for r in rows),
                'last_observation': rows[-1], 'history': rows})
        if time.monotonic()-last_print >= 55:
            try:
                current = load('running')
                row = current.get('last_completed_solve') or {}
                print(f"{now()} {current['status']} step={row.get('step')} nonlinear={row.get('nonlinear_iteration')} linear_solves={current.get('linear_solves_completed')} elapsed_s={current.get('elapsed_s')}", flush=True)
            except (FileNotFoundError, json.JSONDecodeError):
                pass
            last_print = time.monotonic()
        time.sleep(2)

    def script(name):
        with (LOG/(name+'.log')).open('w') as out:
            return subprocess.run([sys.executable, '-B', str(ROOT/'scripts/sv12'/(name+'.py'))],
                                  cwd=ROOT, stdout=out, stderr=subprocess.STDOUT).returncode
    assert script('analyze') == 0, 'Completed run analysis failed; inspect logs/sv1_2/analyze.log'
    reload_code = script('reload_solution')
    if reload_code:
        flow = load('flow_execution')
        flow.update(status='FAIL', reason='RESULT_RELOAD_MISMATCH', accepted_solution_available=False)
        write_json(REPORT/'flow_execution.json', flow)
        (REPORT/'accepted_solution.json').rename(REPORT/'rejected_solution_reload.json')
    assert script('figures') == 0, 'Figure generation failed; inspect logs/sv1_2/figures.log'
    # First audit provides actual evidence for the permanent preservation test.
    script('audit_preservation')
    expected = {
        'tests/test_sv1_mass_balance.py::test_actual_field_mass_conservation',
        'tests/test_sv1_solver_result.py::test_real_solver_converged',
        'tests/test_sv1_solver_result.py::test_required_steady_intervals_reached',
        'tests/test_sv11_short_run_mass.py::test_actual_step_ten_mass_strictly_improves'}
    results = {}
    for label, selected in [('sv12', sorted(str(p.relative_to(ROOT)) for p in (ROOT/'tests').glob('test_sv12_*.py'))), ('all', [])]:
        xml = REPORT/('pytest_'+label+'.xml')
        with (LOG/('pytest_'+label+'.log')).open('w') as out:
            process = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-q', '-p', 'no:cacheprovider',
                                      *selected, '--junitxml='+str(xml)], cwd=ROOT, stdout=out, stderr=subprocess.STDOUT)
        suites = ET.parse(xml).getroot()
        cases = suites.findall('.//testcase')
        failures, errors = [], []
        for case in cases:
            classname = case.attrib['classname']
            name = classname.replace('.', '/')+'.py::'+case.attrib['name']
            if case.find('failure') is not None: failures.append(name)
            if case.find('error') is not None: errors.append(name)
        skipped = sum(case.find('skipped') is not None for case in cases)
        results[label] = {
            'passed': len(cases)-len(failures)-len(errors)-skipped, 'failed': len(failures),
            'skipped': skipped, 'errors': len(errors), 'failures': failures,
            'error_cases': errors, 'process_exit_code': process.returncode,
            'junit_xml': str(xml.relative_to(ROOT))}
    actual_failures = set(results['all']['failures'])
    results.update(historical_expected_failures=sorted(actual_failures & expected),
                   known_historical_failure_ids=sorted(expected),
                   unexpected_regression_failures=sorted(actual_failures-expected),
                   missing_historical_failures=sorted(expected-actual_failures),
                   note='Historical failed-stage tests are retained as real failures, never xfailed or deleted. SV1.2 acceptance is classified separately.')
    write_json(REPORT/'test_results.json', results)
    # Required final audit follows test execution; it must leave old evidence intact.
    script('audit_preservation')
    assert script('write_report') == 0, 'Final report generation failed'
    print((REPORT/'terminal_summary.txt').read_text(), flush=True)
    write_json(REPORT/'completion.json', {'completed_utc': now(), 'status': load('stage_result')['status'],
                                        'report': 'reports/sv1_2/REPORT.md', 'next_stage_started': False})

try:
    main()
except Exception:
    write_json(REPORT/'postprocessing_error.json', {'timestamp': now(), 'traceback': traceback.format_exc()})
    raise
