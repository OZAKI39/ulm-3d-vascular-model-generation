"""Verify RFUS evidence, preserving all historical experiments read-only."""
from pathlib import Path
import sys, json, zipfile, subprocess
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vascular_processing.support_inputs import ROOT, discover, verify, save, digest


def main():
    from vascular_processing.support_rfus_pipeline import preparation
    from vascular_processing.support_rfus_fixtures import generate, actual_straight_test
    from PIL import Image
    import trimesh
    import numpy as np
    c = preparation()
    out = c['out']; validation = out / 'validation'
    s = json.loads((out / 'robust_full_support_summary.json').read_text())
    checks = []
    def check(name, passed, detail=None):
        checks.append(dict(check=name, passed=bool(passed), detail=detail))
    protected = json.loads((validation / 'protected_before.json').read_text())
    check('upstream_and_history_unchanged', verify(protected)['all_unchanged'])
    historical = (validation / 'historical_tests_before.py.txt').read_text()
    current = (ROOT / 'tests/test_s1_8_vascular_support_planner.py').read_text()
    check('historical_test_prefix_unchanged', current.startswith(historical))
    for name in ['robust_underbody_config.json', 'underbody_corridors.csv',
                 'curtain_segments.csv', 'roof_ridges.csv', 'removal_handles.csv',
                 'base_necks.csv', 'underbody_coverage.csv',
                 'support_self_printability.csv', 'segment_removal.csv',
                 'travel_near_miss.csv', 'bambu_validation.json', 'method_comparison.csv',
                 'robust_full_underbody_report.md']:
        p = out / name
        check(name, p.is_file() and p.stat().st_size > 0)
    for name in ['custom_full_underbody_support', 'support_roofs', 'support_curtains',
                 'support_handles', 'support_necks', 'mold_without_support',
                 'mold_with_robust_support']:
        p = out / 'geometry' / (name + '.stl')
        check('geometry_' + name, p.is_file() and p.stat().st_size > 84)
        mesh=trimesh.load_mesh(p)
        check('readable_geometry_' + name,len(mesh.faces)>0 and np.isfinite(mesh.vertices).all())
        check('open_geometry_never_production_' + name,mesh.is_volume or s['status']=='NEEDS_ADJUSTMENT',
              dict(closed_positive_volume=bool(mesh.is_volume)))
    images = sorted((out / 'QC_ROBUST').glob('*.png'))
    check('twenty_QC_images', len(images) == 20)
    for p in images:
        with Image.open(p) as image:
            image.verify()
        check(p.name, p.stat().st_size > 1000)
    expected = {f'RFUS_C{lines}_R{lines}_N{neck}' for lines in (2, 3) for neck in (1, 2)}
    check('all_four_required_initial_candidates', expected <= {r['candidate_id'] for r in s['candidates']})
    for item in s['candidates']:
        folder = Path(item['folder'])
        r = json.loads((folder / 'candidate.json').read_text())
        check(item['candidate_id'] + '_O000_scale_one',
              r['scale_factor'] == 1 and r['orientation']['orientation_id'] == 'O000')
        for label, name in [('full', 'X1C_diagnostic_full_support.3mf'),
                            ('support', 'bambu_support_only/sliced.3mf')]:
            with zipfile.ZipFile(folder / name) as archive:
                settings = json.loads(archive.read('Metadata/project_settings.config'))
                check(item['candidate_id'] + '_' + label + '_CRC', archive.testzip() is None)
                check(item['candidate_id'] + '_' + label + '_actual_Gcode',
                      any(n.endswith('.gcode') for n in archive.namelist()))
                check(item['candidate_id'] + '_' + label + '_auto_off',
                      str(settings.get('enable_support')) in ('0', 'false', 'False'))
                check(item['candidate_id'] + '_' + label + '_X1C_02_ABS',
                      settings.get('printer_model') == 'Bambu Lab X1 Carbon'
                      and [float(x) for x in settings.get('nozzle_diameter', [])] == [.2]
                      and settings.get('filament_type') == ['ABS']
                      and settings.get('print_settings_id') == '0.10mm Standard Vessel'
                      and settings.get('filament_settings_id') == ['Bambu ABS Vessel'])
        u = r['slice']['underbody']
        check(item['candidate_id'] + '_complete_denominator', u['sample_count'] == 1990)
        check(item['candidate_id'] + '_no_native_support', r['slice']['native_support_toolpaths'] == 0)
        check(item['candidate_id'] + '_honest_status',
              (r['status'] == 'NEEDS_ADJUSTMENT') == bool(r['failed_gates']))
    check('no_false_production', s['production_selected'] is None
          if s['status'] == 'NEEDS_ADJUSTMENT' else s['production_selected'] is not None)
    check('no_printer_command', s['printer_commands_sent'] is False)
    fixtures = generate(ROOT / 'tests/data/support/rfus', c['profile'])
    check('fifteen_closed_fixtures', len(fixtures) >= 15 and all(r['closed'] for r in fixtures))
    fixture_result = validation / 'actual_straight_tube/result.json'
    actual = json.loads(fixture_result.read_text()) if fixture_result.is_file() else actual_straight_test(c)
    if not actual.get('real_sparse_control'):actual=actual_straight_test(c)
    with zipfile.ZipFile(actual['archive']) as archive:
        check('actual_straight_fixture_archive_CRC', archive.testzip() is None)
        check('actual_straight_fixture_contains_Gcode', any(n.endswith('.gcode') for n in archive.namelist()))
    with zipfile.ZipFile(actual['sparse_archive']) as archive:
        check('actual_sparse_control_archive_CRC', archive.testzip() is None)
        check('actual_sparse_control_contains_Gcode', any(n.endswith('.gcode') for n in archive.namelist()))
    check('actual_long_horizontal_full_underbody_and_sparse_negative_control', actual['passed'],
          {k: v for k, v in actual.items() if k != 'full_underbody'})
    command = [str(ROOT / 'run_python.sh'), '-m', 'pytest', '-q',
               'tests/test_s1_8_vascular_support_planner.py',
               '--junitxml=' + str(validation / 'pytest.xml')]
    result = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    (validation / 'pytest.log').write_text(result.stdout)
    check('full_persistent_regression', result.returncode == 0, result.stdout.splitlines()[-1:])
    protection = verify(protected)
    check('final_source_protection', protection['all_unchanged'])
    report = dict(artifact_workflow_passed=all(x['passed'] for x in checks),
                  production_status=s['status'], checks=checks, source_protection=protection,
                  code_sha256={str(p.relative_to(ROOT)): digest(p)
                               for p in (ROOT / 'vascular_processing').glob('support_rfus*.py')})
    save(validation / 'artifact_verification.json', report)
    marker = '\n## 永久测试与独立切片验证\n'
    path = out / 'robust_full_underbody_report.md'
    text = path.read_text().split(marker)[0]
    text += marker + '\n' + (
        f"独立 30 mm 水平圆管已执行真实 Bambu 切片，下托覆盖 "
        f"{actual['full_underbody']['percent']:.6f}%；稀疏负对照的最长缺口 "
        f"{actual['sparse_negative_control']['max_uncovered_underbody_arc_mm']:.6f} mm。"
        "该试件仅验证真实刀路下托指标，不能替代完整模型的拆除或物理验证。\n\n"
        f"永久回归结果：{result.stdout.splitlines()[-1]}。"
        f"工件与证据完整性检查通过：{report['artifact_workflow_passed']}；"
        f"完整模型状态仍为 {s['status']}。\n")
    path.write_text(text)
    print(json.dumps({k: v for k, v in report.items() if k not in ('checks', 'code_sha256')},
                     ensure_ascii=False, indent=2))
    if not report['artifact_workflow_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
