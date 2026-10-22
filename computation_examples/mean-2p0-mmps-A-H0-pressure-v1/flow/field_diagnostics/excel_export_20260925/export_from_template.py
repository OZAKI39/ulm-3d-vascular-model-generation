"""Export the current H0 residuals into the supplied workbook without restyling.

Only numeric worksheet data, native chart caches, provenance text and document
metadata are replaced. The source workbook and all CFD evidence stay read-only.
Run with the project Python environment (lxml required).
"""
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile
import csv
import hashlib
import json
import math
import re
import sys
from lxml import etree as E

HERE = Path(__file__).resolve().parent
CASE = HERE.parents[1]
ROOT = CASE.parents[1]
TEMPLATE = Path('/home/lzy/projects/A_Global_Transient_Evolution.xlsx')
DEST = TEMPLATE.with_name('A_Global_Transient_Evolution_A_H0.xlsx')
sys.path.insert(0, str(ROOT / 'scripts/sv13q'))
from flow_parser import parse_solver_log, ROW

NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'c': 'http://schemas.openxmlformats.org/drawingml/2006/chart'}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def write_csv(path, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def xml_bytes(root):
    return E.tostring(root, encoding='UTF-8', xml_declaration=True, standalone=True)


def main():
    assert not DEST.exists(), 'Never overwrite a prior workbook'
    inputs = [TEMPLATE, CASE/'run/solver.log', CASE/'run/solver.xml',
              CASE/'run/1-procs/histor.dat', CASE/'reports/execution.json',
              CASE/'reports/physics_validation_H0.json', CASE/'frozen_flow/manifest.json',
              ROOT/'scripts/sv13q/flow_parser.py',
              ROOT/'vendor/svMultiPhysics_stage_q/Code/Source/solver/Integrator.cpp',
              ROOT/'vendor/svMultiPhysics_stage_q/Code/Source/solver/petsc_impl.cpp']
    hashes = {str(p): sha(p) for p in inputs}
    config = E.parse(str(CASE/'run/solver.xml'))
    equation = config.find('.//Add_equation[@type="fluid"]')
    dt = float(config.findtext('.//Time_step_size'))
    tolerance = float(equation.findtext('Tolerance'))
    minimum = int(equation.findtext('Min_iterations'))
    maximum = int(equation.findtext('Max_iterations'))
    assert config.findtext('.//Continue_previous_simulation').lower() == 'false'
    boundary = {n: float(equation.findtext(f"Add_BC[@name='{n}']/Value"))
                for n in ['OUTLET_01', 'OUTLET_02', 'OUTLET_03']}
    execution = json.loads((CASE/'reports/execution.json').read_text())
    physics = json.loads((CASE/'reports/physics_validation_H0.json').read_text())
    frozen = json.loads((CASE/'frozen_flow/manifest.json').read_text())
    assert execution['status'] == physics['status'] == frozen['status'] == 'PASS'
    assert execution['log_sha256'] == sha(CASE/'run/solver.log')
    assert frozen['input_configuration_sha256'] == sha(CASE/'run/solver.xml')
    flow_name = 'steady_flow_mean_2p0_mmps_A_H0.vtu'
    assert frozen['files'][flow_name]['sha256'] == '064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4'
    log = (CASE/'run/solver.log').read_text()
    history = parse_solver_log(log, dt)
    assert history == execution['history']
    assert not any(history[k] for k in ['unparsed_rows', 'unassigned_warnings',
        'unassigned_petsc_monitor', 'failed_linear_solves', 'nonlinear_failure_messages', 'recovered_attempts'])
    events = [dict(re.findall(r'(\w+)=([^\s]+)', s)) for s in log.splitlines()
              if s.startswith('SV13Q_BEGIN ')]
    assert all(e['attempt'] == '0' for e in events)
    # Independent NS table comparison with the separately saved native history.
    def native_rows(text):
        return [m.groupdict() for line in text.splitlines() if (m := ROW.match(line))]
    assert native_rows(log) == native_rows((CASE/'run/1-procs/histor.dat').read_text())
    solves = history['linear_solves']
    assert len(solves) == len(events)
    groups = {}
    for s in solves:
        groups.setdefault(s['step'], []).append(s)
        assert s['linear_converged'] and s['residual_trustworthy']
        assert s['petsc_reason']['reason'].startswith('CONVERGED_')
        assert s['petsc_monitor'][0]['iteration'] == 0
        assert all(m['finite'] for m in s['petsc_monitor'])
    assert list(groups) == list(range(1, execution['final_step'] + 1))
    reference = solves[0]['petsc_monitor'][0]['residual_norm']
    assert solves[0]['step'] == solves[0]['nonlinear_iteration'] == 1
    iteration_rows, step_rows, transitions = [], [], []
    for n, group in groups.items():
        assert [s['nonlinear_iteration'] for s in group] == list(range(1, len(group)+1))
        start = group[0]['petsc_monitor'][0]['residual_norm']
        for k, s in enumerate(group, 1):
            monitor = s['petsc_monitor'][0]
            norm = monitor['residual_norm']
            assert norm > 10 * sys.float_info.epsilon**2
            assert math.isclose(norm, monitor['true_residual_norm'], rel_tol=2e-12)
            global_ratio, local_ratio = norm/reference, norm/start
            assert math.isclose(global_ratio, s['nonlinear_Ri_over_R0'], rel_tol=5.1e-4)
            assert math.isclose(local_ratio, s['nonlinear_Ri_over_R1'], rel_tol=5.1e-4)
            stop = k >= maximum or (k >= minimum and
                                    (global_ratio <= tolerance or local_ratio <= tolerance))
            assert stop == (k == len(group)) and k < maximum
            event = events[s['linear_solve_index']-1]
            assert int(event['step']) == n and int(event['logical']) == s['linear_solve_index']
            iteration_rows.append(dict(step=n, iteration=k, initial_KSP_norm=norm,
                global_ratio=global_ratio, step_relative_ratio=local_ratio,
                printed_Ri_over_R0=s['nonlinear_Ri_over_R0'], printed_Ri_over_R1=s['nonlinear_Ri_over_R1'],
                initial_KSP_log_line=monitor['line'], NS_log_line=s['log_line'],
                stopping_rule_satisfied=stop))
        end = group[-1]['petsc_monitor'][0]['residual_norm']
        step_rows.append(dict(step=n, time_s=n*dt, iterations=len(group), start_norm=start,
            end_norm=end, start_global=start/reference, end_global=end/reference,
            end_over_start=end/start, initial_KSP_log_line=group[0]['petsc_monitor'][0]['line']))
        if n > 1 and len(group) != len(groups[n-1]):
            transitions.append(dict(step=n, previous_iterations=len(groups[n-1]), current_iterations=len(group)))
    plot_rows = [dict(Step=r['step'], Start_global=r['start_global'],
                      Global_tolerance=tolerance, Time_s=r['time_s']) for r in step_rows]
    write_csv(HERE/'Plot_Data.csv', plot_rows)
    write_csv(HERE/'nonlinear_step_audit.csv', step_rows)
    write_csv(HERE/'nonlinear_iteration_audit.csv', iteration_rows)
    audit = dict(case=str(CASE), all_pass=True, source_hashes=hashes,
        flow_sha256=frozen['files'][flow_name]['sha256'], outlet_BC_Pa=boundary,
        definition='Start_global = first logged initial KSP norm of this time step / first logged initial KSP norm of step 1; fixed reference throughout run',
        reference_norm=reference, flow_solver_dt_s=dt, time_steps=len(step_rows), nonlinear_records=len(solves),
        nonlinear_tolerance=tolerance, minimum_iterations=minimum, maximum_iterations=maximum,
        stopping_rule='k >= max OR (k >= min AND (R/R_ref <= tol OR R/R_step_start <= tol))',
        iteration_transitions=transitions, maximum_iteration_cap_never_reached=True,
        all_steps_stop_at_first_eligible_record=True, native_histor_matches_solver_log=True,
        parsed_history_matches_existing_execution_audit=True, numerical_floor_established=False,
        raw_monitor_precision_preserved=True, time_axis_is_flow_solver_time_not_particle_time=True,
        first_row=plot_rows[0], last_row=plot_rows[-1])
    write_json(HERE/'NONLINEAR_RESIDUAL_AUDIT.json', audit)

    with ZipFile(TEMPLATE) as archive:
        infos = archive.infolist()
        original = {i.filename: archive.read(i.filename) for i in infos}
    parts = dict(original)
    strings_root = E.fromstring(parts['xl/sharedStrings.xml'])
    strings = [''.join(s.itertext()) for s in strings_root]
    sheet = E.fromstring(parts['xl/worksheets/sheet1.xml'])
    cells = {c.get('r'): c for c in sheet.findall('.//s:sheetData/s:row/s:c', NS)}
    assert [strings[int(cells[f'{c}1'].find('s:v', NS).text)] for c in 'ABCD'] == list(plot_rows[0])
    assert len(sheet.findall('s:sheetData/s:row', NS)) == len(plot_rows)+1
    for i, row in enumerate(plot_rows, 2):
        for col, value in zip('ABCD', row.values()):
            cell = cells[f'{col}{i}']
            assert cell.get('t', 'n') == 'n'
            cell.find('s:v', NS).text = repr(value)
    parts['xl/worksheets/sheet1.xml'] = xml_bytes(sheet)

    notes = E.fromstring(parts['xl/worksheets/sheet3.xml'])
    note_cells = {c.get('r'): c for c in notes.findall('.//s:sheetData/s:row/s:c', NS)}
    all_string_refs = Counter(int(c.find('s:v', NS).text)
        for n, raw in original.items() if re.fullmatch(r'xl/worksheets/sheet\d+.xml', n)
        for c in E.fromstring(raw).findall('.//s:c', NS) if c.get('t') == 's')
    def note(address, value):
        index = int(note_cells[address].find('s:v', NS).text)
        assert all_string_refs[index] == 1, 'Do not change strings shared with another cell'
        item = strings_root[index]
        assert len(item) == 1 and item[0].tag == '{'+NS['s']+'}t'
        item[0].text = value
    note('B2', f"新边界条件流场：{CASE.name}；O1/O2/O3 = " +
         '/'.join(f'{x:.6f}' for x in boundary.values()) + ' Pa；对应残差图 A。')
    note('B8', f'R_ref = 第 1 时间步第 1 次非线性记录的范数 = {reference:.13e}；全程固定，不逐步重新归一化。')
    note('A11', '非线性迭代变化')
    note('B11', '；'.join(f"第 {t['step']} 步：{t['previous_iterations']} → {t['current_iterations']} 次迭代" for t in transitions)
         + '。全部时间步满足停止条件，未触及迭代上限。')
    note('B15', f'对应流场求解时间，单位 s；Time_s = Step × {dt:.16e}。这里使用 FEM 求解时间步，不使用微泡轨迹的 1.0 ms 时间步。')
    sources = [HERE/'nonlinear_step_audit.csv', CASE/'run/solver.log', CASE/'run/solver.xml',
               HERE/'NONLINEAR_RESIDUAL_AUDIT.json', CASE/'reports/execution.json']
    for i, p in enumerate(sources, 16):
        note(f'A{i}', '来源 SHA256：'+p.name)
        note(f'B{i}', sha(p))
    parts['xl/sharedStrings.xml'] = xml_bytes(strings_root)

    chart = E.fromstring(parts['xl/charts/chart1.xml'])
    for ref in chart.findall('.//c:numRef', NS):
        formula = ref.findtext('c:f', namespaces=NS)
        m = re.fullmatch(r'Plot_Data!\$([A-D])\$2:\$\1\$(\d+)', formula)
        assert m and int(m[2]) == len(plot_rows)+1
        column = 'ABCD'.index(m[1])
        cache = ref.find('c:numCache', NS)
        points = cache.findall('c:pt', NS)
        assert len(points) == int(cache.find('c:ptCount', NS).get('val')) == len(plot_rows)
        for i, point in enumerate(points):
            assert int(point.get('idx')) == i
            point.find('c:v', NS).text = repr(list(plot_rows[i].values())[column])
    parts['xl/charts/chart1.xml'] = xml_bytes(chart)
    core = E.fromstring(parts['docProps/core.xml'])
    core.find('{http://purl.org/dc/terms/}modified').text = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    core.find('{http://purl.org/dc/elements/1.1/}subject').text = 'New H0 outlet-pressure case; supplied template styling preserved.'
    parts['docProps/core.xml'] = xml_bytes(core)
    with ZipFile(DEST, 'x') as archive:
        for info in infos:
            archive.writestr(info, parts[info.filename])

    # Reopen the actual deliverable; validate every cell and every chart cache.
    with ZipFile(DEST) as archive:
        assert archive.testzip() is None
        saved = {n: archive.read(n) for n in archive.namelist()}
    assert saved.keys() == original.keys()
    allowed = {'xl/worksheets/sheet1.xml', 'xl/sharedStrings.xml', 'xl/charts/chart1.xml', 'docProps/core.xml'}
    changed = {n for n in saved if saved[n] != original[n]}
    assert changed <= allowed
    numeric_checks = 0
    saved_cells = {c.get('r'): c for c in E.fromstring(saved['xl/worksheets/sheet1.xml']).findall('.//s:sheetData/s:row/s:c', NS)}
    for i, row in enumerate(plot_rows, 2):
        for col, expected in zip('ABCD', row.values()):
            c = saved_cells[f'{col}{i}']
            assert c.get('t', 'n') == 'n'
            assert float(c.findtext('s:v', namespaces=NS)) == expected
            numeric_checks += 1
    chart_check = E.fromstring(saved['xl/charts/chart1.xml'])
    for ref in chart_check.findall('.//c:numRef', NS):
        col = 'ABCD'.index(re.search(r'!\$([A-D])', ref.findtext('c:f', namespaces=NS))[1])
        for point in ref.findall('c:numCache/c:pt', NS):
            assert float(point.findtext('c:v', namespaces=NS)) == list(plot_rows[int(point.get('idx'))].values())[col]
    def without_cache(raw):
        root = E.fromstring(raw)
        for e in root.findall('.//c:numCache', NS): e.getparent().remove(e)
        return E.tostring(root, method='c14n')
    assert without_cache(saved['xl/charts/chart1.xml']) == without_cache(original['xl/charts/chart1.xml'])
    def without_numeric_values(raw):
        root = E.fromstring(raw)
        for c in root.findall('.//s:sheetData/s:row/s:c', NS):
            if c.get('t', 'n') == 'n': c.remove(c.find('s:v', NS))
        return E.tostring(root, method='c14n')
    assert without_numeric_values(saved['xl/worksheets/sheet1.xml']) == without_numeric_values(original['xl/worksheets/sheet1.xml'])
    assert all(sha(p) == digest for p, digest in hashes.items())
    verification = dict(status='PASS', output=str(DEST), output_sha256=sha(DEST),
        template=str(TEMPLATE), template_sha256=hashes[str(TEMPLATE)], template_unchanged=True,
        data_rows=len(plot_rows), numeric_cells_verified=numeric_checks,
        chart_cache_values_verified=True, chart_style_and_formulas_unchanged=True,
        worksheet_layout_and_cell_styles_unchanged=True, unchanged_zip_parts=len(saved)-len(changed),
        changed_zip_parts=sorted(changed), source_files_unchanged=True,
        reference_norm=reference, first_Start_global=plot_rows[0]['Start_global'],
        last_Start_global=plot_rows[-1]['Start_global'], flow_solver_dt_s=dt,
        nonlinear_iteration_transitions=transitions,
        provenance_directory=str(HERE), script_sha256=sha(Path(__file__)))
    write_json(HERE/'EXPORT_VALIDATION.json', verification)
    print(json.dumps(verification, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
