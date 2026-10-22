"""Small acceptance checks for the current mouse workflow; no full simulation.

Run using this project's editable Python environment from any working directory.
Existing tests and scientific functions are reused without changing tolerances.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import traceback
import xml.etree.ElementTree as ET

from mouse_vessel_flow.paths import PROJECT_ROOT, new_output
from mouse_vessel_flow.provenance import verify_source_integrity, protected_originals
from mouse_vessel_flow.network import validate_design
from mouse_vessel_flow.validation import field_check, wss_check
from mouse_vessel_flow.native_environment import native_status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, help='New output directory; never overwrite existing evidence')
    args = parser.parse_args()
    out = new_output(args.output)
    checks = {}

    def check(name, function):
        try:
            checks[name] = function()
        except Exception as exc:
            checks[name] = {'status': 'FAIL', 'error': repr(exc), 'traceback': traceback.format_exc()}
        print(name + ': ' + checks[name]['status'], flush=True)

    check('current_source', lambda: {'status': 'PASS', 'count': len(verify_source_integrity())})
    check('original_source_evidence', lambda: {'status': 'PASS', 'count': len(protected_originals())})
    check('small_pressure_network', validate_design)
    check('current_saved_field', lambda: field_check('roi_only'))
    check('current_saved_wall_quantity', wss_check)
    check('installed_native_CFD', lambda: native_status('native_cfd'))
    check('installed_native_mesher', lambda: native_status('native_mesher'))
    check('installed_surface_tools', lambda: native_status('vmtk'))
    check('installed_lumen_tool', lambda: native_status('ultraliser'))

    # Select the current case's source/entry test explicitly; keep its assertions.
    base = 'tests/test_refactor_equivalence.py::'
    nodes = ['tests/network/test_roi_only.py'] + [base + name for name in (
        'test_full_small_scene_outputs_match_original_bytes',
        'test_native_handoff_kernel_matches_original_frozen_wall_bytes',
        'test_retained_geometry_is_the_flow_roi',
        'test_frozen_output_guard_rejects_write',
        'test_embedded_tetgen_keeps_original_body_and_legacy_syntax',
        'test_real_source_admission_and_job_binding[roi_only]',
    )]
    command = [sys.executable, '-m', 'pytest', '-q', *nodes, '--junitxml=' + str(out / 'pytest.xml')]
    with (out / 'pytest.log').open('x') as stream:
        run = subprocess.run(command, cwd=PROJECT_ROOT, stdout=stream, stderr=subprocess.STDOUT)
    cases = list(ET.parse(out / 'pytest.xml').getroot().iter('testcase')) if (out / 'pytest.xml').exists() else []
    failed = sum(c.find('failure') is not None or c.find('error') is not None for c in cases)
    skipped = sum(c.find('skipped') is not None for c in cases)
    checks['small_tests'] = {
        'status': 'PASS' if run.returncode == 0 and cases and not skipped else 'FAIL',
        'passed': len(cases) - failed - skipped, 'failed': failed, 'skipped': skipped,
        'returncode': run.returncode, 'command': command,
    }
    print('small_tests: ' + checks['small_tests']['status'], flush=True)
    report = {
        'current_workflow': 'roi_only',
        'tested_components_status': 'PASS' if all(c['status'] == 'PASS' for c in checks.values()) else 'FAIL',
        'checks': checks,
        'native_mesher_execution': checks['installed_native_mesher'],
        'native_CFD_execution': checks['installed_native_CFD'],
        'full_simulation': {'status': 'NOT_REQUIRED_BY_USER', 'executed': False},
        'current_test_inventory': 'validation/TEST_STATUS.json',
        'meaning': 'Small tests accept the components they exercise. Native execution evidence is separately recorded and hash verified; no full mouse run is required.',
    }
    (out / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    return 0 if report['tested_components_status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
