#!/usr/bin/env python3
"""Derive sample-level ROI cuts from frozen predictions, then run native VMD."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import networkx as nx
import numpy as np

from vascular_processing.boundary_review_export import load_cache, protect_inputs
from vascular_processing.mevo_refinement import (
    RULES, branch_edges, export_subset, frozen_upstream, graph_length, prune_diameter,
    refine_component, statistics, table)
from vascular_processing.mevo_refinement_figures import export_figures
from vascular_processing.swc_export import read_source, validate_tree
from vascular_processing.topbrain_qc import sha256


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write('\n')


def protection(cache, output):
    snapshot = dict(cache.snapshot)
    previous = cache.directory.parent.parent / 'protection_before.json'
    for relative, digest in json.loads(previous.read_text()).items():
        path = ROOT / relative
        if not path.is_file() or sha256(path) != digest:
            raise ValueError('Pre-existing source protection mismatch: ' + relative)
        snapshot[str(path)] = digest
    # Includes existing review ZIP, all 295 review files and all old models.
    for path in cache.directory.rglob('*'):
        if path.is_file() and not path.is_relative_to(output):
            snapshot[str(path)] = sha256(path)
    return snapshot


def model_one(path, output):
    import pyvista as pv
    from vascular_processing.pipeline import Options, process_file
    started = time.monotonic()
    print('VASCULARMD START', path.name, flush=True)
    report = process_file(path, output, Options())
    failures = report['failed_components']
    called = not any(f['stage'] in {'input_validation', 'native_load_and_preprocessing'} for f in failures)
    model_success = called and not any(f['stage'] == 'native_model' for f in failures)
    result = dict(status=report['status'], model_network_called=called, model_network_success=model_success,
        native_model_call=report['native_model_call'], topology_preserved=report.get('topology_preserved'),
        outputs=report['outputs'], failures=[{k: f[k] for k in ('stage', 'error_type', 'message')} for f in failures],
        warnings=report['warnings'], elapsed_seconds=time.monotonic() - started,
        qc_path=str(output / (path.stem + '_vmd_qc.json')))
    vtk = result['outputs'].get('surface_vtk')
    if vtk:
        mesh = pv.read(vtk).extract_surface(algorithm='dataset_surface').triangulate()
        stl = output / (path.stem + '_surface.stl')
        mesh.save(stl)
        check = pv.read(stl)
        if not check.n_cells or not np.isfinite(check.points).all():
            raise ValueError('Invalid STL read-back')
        result['outputs']['surface_stl'] = str(stl)
        result['stl_qc'] = dict(points=check.n_points, cells=check.n_cells, open_edges=check.n_open_edges,
            diagnostic_only=report.get('surface', {}).get('diagnostic_only', False),
            manufacturing_validated=False)
    print('VASCULARMD END', path.name, result['status'], result['failures'], flush=True)
    return result


def report_text(summary):
    totals = summary['totals']
    lines = ['# 冻结 NN 结果支持的 MeVO 近端边界精修与持续性管径筛选', '',
        '## 摘要', '',
        '本研究读取 BG001 右侧大脑中动脉的既有多供体 NN 预测，将整分支边界细化到原始 SWC 采样点，并对精修结果实施远端持续细径筛选。配准、供体排序、NN 查询、UNKNOWN 校准和分支概率聚合均未重新运行。两层 SWC 分别保留，以区分语义证据与研究设定的管径条件。', '',
        f"原候选总长度为 {totals['original_length_mm']:.6f} mm，语义精修后为 {totals['semantic_length_mm']:.6f} mm，直径筛选后为 {totals['diameter_length_mm']:.6f} mm。直径阶段删除 {totals['diameter_removed_nodes']} 个原始节点、{totals['diameter_removed_length_mm']:.6f} mm 中心线，确认 {totals['sustained_cuts']} 处持续细径截断。", '',
        '## 方法', '',
        '所有坐标、半径和 TYPE 均取自冻结的原始 BG001.CNG.swc。SWC 第六列是半径，故最低直径 0.75 mm 对应半径 0.375 mm。派生 SWC 仅删除前缀或子树、重新编号并设定新根；未改变原始几何或数值。VascularMD 的平滑结果另存，不反向参与筛选。', '',
        '近端候选必须具有 pMeVO≥0.60。在沿原始中心线弧长的闭区间 [起点, 起点+5 mm] 内，至少三个有效原始采样点中不少于 80% 满足相同条件，且窗口中不得出现 pM1≥0.80。分支不足 5 mm 时使用剩余长度。有效点要求冻结标签不为 UNKNOWN、有效供体数至少为 3，并具有有限的概率、支持距离和一致性。若候选之后至该拓扑分支末端仍出现跨越至少 3 mm 的连续 pM1≥0.60 采样序列，则拒绝该候选。反转检查也纳入 UNKNOWN 点的有限 M1 概率，不让低支持标签掩盖反转；仅非有限 M1 概率打断该序列。稀疏采样不被补点。', '',
        'part03 根据本轮指定范围额外搜索 UNKNOWN 分支 94_112；这不会改写原分支标签。每个女儿入口独立检验。相同的保留节点可合并；不相连的入口分别导出单根 SWC。表中位移沿原始流向为正，向上游扩展为负。', '',
        '直径判据使用最短的、覆盖至少 3 mm 的原始样本窗口，末端取首次到达或越过 3 mm 的原始点；实际跨度逐条记录，因此稀疏采样时可以大于 3 mm。窗口首点必须低于 0.75 mm，至少有两个低径样本且低径样本比例不低于 80%，避免单个低径点触发删除。此处“持续”按物理跨度和样本占比定义，允许最多 20% 样本波动，不表示对连续管腔作了测量。若随后严格 5 mm 范围内出现连续原始高径样本跨越至少 3 mm，则记录 DIAMETER_REBOUND 并否决该候选。', '',
        '逐条根至末端路径评估，并按距精修根的距离从近到远裁切。共享候选任一路径具有回弹证据即否决；确认截断时保留截断点，删除全部后代，因此不会留下远端孤岛。已被近端截断覆盖的候选仅保留审计记录，不重复计数。分支起始接点可作为上游末端保留，其后整段及子树删除；若精修根本身低于阈值，则整个派生组件不进入直径版。', '',
        '**直径筛选是 MeVO 研究及制造目标的操作性筛选，不是解剖学 M3 终点，也不是 M3→M4 标注。** 本规则有意容忍孤立小半径点，并保留低于阈值的截断端点，因此最终文件不保证每一个点都达到 0.75 mm；不能把该文件名理解为逐点硬阈值。', '',
        '## 结果', '',
        '| 组件 / 入口 | 原始根 | 精修根 | 位移 mm | 状态 |', '|---|---:|---:|---:|---|']
    for r in summary['proximal_refinement']:
        shift = r['distance_from_original_boundary_mm']
        lines.append(f"| part{r['component']:02d} / {r['refined_cut_branch']} | {r['original_boundary_node']} | {r['refined_cut_original_swc_id']} | {shift:.6f} | {r['status']} |" if shift is not None else f"| part{r['component']:02d} / {r['refined_cut_branch']} | {r['original_boundary_node']} | — | — | {r['status']} |")
    lines += ['', '| 组件 | 原长度 mm | 语义长度 mm | 管径筛选后 mm | 删除节点 | 删除旧末端 | 持续截断 |', '|---|---:|---:|---:|---:|---:|---:|']
    for c in summary['components']:
        lines.append(f"| part{c['component']:02d} | {c['original']['length_mm']:.6f} | {c['semantic']['length_mm']:.6f} | {c['diameter']['length_mm']:.6f} | {c['removed_node_count']} | {c['removed_original_terminal_count']} | {c['sustained_cuts']} |")
    lines += ['', '“删除旧末端”计数的是语义版原有末端中被删除的数量，不是前后末端数相减；截断会产生新末端。完全删除分支依据该原始分支所有保留边均消失判定；部分删除另外计数。分支等效数按各原始分支保留长度/原始完整长度求和，不能与连通组件数混用。', '',
        '| 组件 | 语义 D 最小/中位/最大 mm | 筛选后 D 最小/中位/最大 mm | 保留分支等效数 | 完全/部分删除分支 |', '|---|---|---|---:|---|']
    for c in summary['components']:
        fmt = lambda s: '/'.join(f'{s[k]:.4f}' if s[k] is not None else '—' for k in ('min_diameter_mm','median_diameter_mm','max_diameter_mm'))
        lines.append(f"| part{c['component']:02d} | {fmt(c['semantic'])} | {fmt(c['diameter'])} | {c['diameter']['branch_equivalent_count']:.6f} | {len(c['fully_removed_branches'])}/{len(c['partly_removed_branches'])} |")
    low_samples = sum(c['diameter']['retained_subthreshold_samples'] for c in summary['components'])
    lines += ['', f"直径阶段完全删除 {totals['fully_removed_branches']} 条原始分支，部分裁短 {totals['partly_removed_branches']} 条，删除 {totals['removed_original_terminals']} 个语义版旧末端。最终仍含 {low_samples} 个 D<0.75 mm 原始点，包含低径截点及未形成持续细径窗的波动样本。真实数据中的回弹否决为 {totals['rebound_vetoes']} 次；合成测试另外验证回弹条件确实能阻止截断。"]
    lines += ['', '低语义支持分支仅另列 SEMANTIC_BORDERLINE 标记，未据此删除：' + '，'.join(
        f"{b['branch_id']}（mean pMeVO={b['p_MeVO']:.6f}）" for b in summary['borderline_branches']) + '。', '',
        '完全被直径规则排除的派生组件：' + '，'.join(summary['excluded_derived_components']) + '。', '',
        '## 原生建模与独立验证', '',
        '仅在全部派生 SWC 通过单根、无环、父节点有效、数值有限且半径为正，以及逐点几何/TYPE 和逐边方向的读回验证后，调用既有 VascularMD 网络流程。使用 radius_model=True、criterion="AIC"，保持现有严格拓扑检查；不为通过建模而修改半径或阈值。', '',
        '| 输入 | 自然单入口 | model_network 调用/完成 | 最终状态 | VTK/STL |', '|---|---|---|---|---|']
    for name, m in summary['models'].items():
        lines.append(f"| {name} | {m.get('single_inlet_stem')} | {m['model_network_called']}/{m['model_network_success']} | {m['status']} | {'有' if m['outputs'].get('surface_vtk') else '无'}/{'有' if m['outputs'].get('surface_stl') else '无'} |")
        for f in m['failures']:
            lines.append(f"\n{name} 的失败发生于 `{f['stage']}`：{f['message']}。语义组件照常保留。\n")
    if 'roi_part03_semantic_refined' in summary['models']:
        sem = summary['models']['roi_part03_semantic_refined']
        dia = summary['models'].get('roi_part03_d075', {})
        lines += ['', f"重点组件 part03 从原始节点 1697 出发，保留 94_112 的真实后缀，形成通向 1702 分叉点的自然单入口；未添加 M1 context 或虚拟血管。语义版原生建模状态为 {sem['status']}，直径版为 {dia.get('status', '未运行')}。"]
    lines += ['', '对于严格语义版根节点仍分叉的输入，先保留实际失败记录，再仅以真实 BraVa 上游分支另存建模输入及 context 映射；该段不计入语义或直径 ROI。所有 VTK/STL 均是建模结果或明确标注的诊断结果，未宣称已通过制造或 CFD 验证。', '',
        '## 可复核性与局限', '',
        '原始预测、五个候选 ROI、既有建模输出、外部审查包及 VascularMD 上游文件均受哈希快照保护。全部源哈希见 protection_before.json，结束验证见 protection_verification.json。代码运行时禁用上游算法入口，未重新配准或查询 NN。', '',
        '固定规则依赖既有 NN 语义与原始半径质量；它提供可审计的计算边界，不替代人工解剖标注。有限采样使物理窗端点离散化，短段不足三个有效点时不会被强行判定。对近端小于阈值的整段排除来自本轮明确规则，应与持续远端截断分别理解。', '',
        '输出中的 proximal_refinement.csv 给出每个入口证据；diameter_pruned_segments.csv 给出截断、整段排除及回弹否决；decision_audit.json 保存所有候选判断。mappings/ 同时保留节点和有向边映射。QC/ 包含五组近端/管径对比及三张同相机全局图；图中文字明确标出截点实际 D 和 r。', '',
        '测试结果另见 QC/test_results.json、对应日志和 JUnit XML。', '']
    if summary.get('validation'):
        v = summary['validation']
        lines += ['## 实际测试与数据完整性', '',
            f"完整回归 {v['full_suite_passed']} 项通过（既有 {v['stable_passed']} 项、新增 {v['new_passed']} 项），零最终测试失败；最后的规则复核后，新增测试再次全部通过。原有 TopBrain 测试产生 {v['stable_warnings']} 条 scikit-image/NumPy 弃用提示。原生建模中的分叉根拒绝属于 part05 严格输入的真实失败，已与测试失败分开记录。", '',
            f"初次 STL 转换记录 {v['initial_stl_future_warnings']} 条 PyVista 将来默认参数变化提示。现已显式固定为本次实际使用的 dataset_surface 算法，并对九个 VTK/STL 的读回重新验证，未改写任何模型文件。", '',
            f"结束时 {v['protected_files']} 个受保护文件哈希全部一致；9 个派生 ROI 通过最终规则重放与节点、边的逐项核对。下列关键源 SHA256 可独立核验：", '',
            '| 源文件 | SHA256 |', '|---|---|']
        lines += [f'| {name} | `{digest}` |' for name,digest in v['key_hashes'].items()]
        lines += ['']
    return '\n'.join(lines)


def verify_results(input_dir, output):
    """Read-only replay of cut decisions and exact exported node/edge mapping."""
    from vascular_processing.boundary_review_export import read_csv
    cache = load_cache(input_dir)
    summary = json.loads((output/'refinement_summary.json').read_text())
    snapshot = json.loads((output/'protection_before.json').read_text())
    source = read_source(Path(cache.manifest['source']['raw_source']))
    edges = branch_edges(cache)
    checked = []
    proximal_rows, diameter_rows = [], []
    with protect_inputs(snapshot), frozen_upstream():
        for component in cache.components:
            _, graphs, rows, _ = refine_component(cache, component, source)
            proximal_rows.extend(rows)
            for index, graph in enumerate(graphs, 1):
                derived = f'part{component:02d}' + (f'_{index:02d}' if len(graphs)>1 else '')
                pruned, rows, _ = prune_diameter(graph, cache.coordinates, edges, component, derived)
                diameter_rows.extend(rows)
                for level, expected in [('semantic_refined', graph), ('diameter_075', pruned)]:
                    exports = [e for e in summary['exports'] if e['derived_component']==derived and e['level']==level]
                    if not expected:
                        assert not exports, 'Excluded component unexpectedly has a filtered SWC'
                        continue
                    assert len(exports)==1
                    exported = exports[0]
                    path = Path(exported['path'])
                    assert sha256(path)==exported['sha256'], 'Derived SWC hash changed'
                    actual = read_source(path).graph
                    rows = read_csv(output/'mappings'/(path.stem+'_nodes.csv'))
                    mapping = {int(r['roi_node_id']):int(r['original_swc_id']) for r in rows}
                    assert set(mapping)==set(actual) and set(mapping.values())==set(expected)
                    assert len(set(mapping.values()))==len(mapping)
                    assert {(mapping[a],mapping[b]) for a,b in actual.edges}==set(expected.edges)
                    for n,original in mapping.items():
                        np.testing.assert_array_equal(actual.nodes[n]['coords'],source.graph.nodes[original]['coords'])
                        assert actual.nodes[n]['swc_type']==source.graph.nodes[original]['swc_type']
                    checked.append(str(path))
        assert proximal_rows==summary['proximal_refinement'], 'Proximal decisions differ'
        assert diameter_rows==summary['diameter_pruning'], 'Diameter decisions differ'
    changed = [p for p,h in snapshot.items() if not Path(p).is_file() or sha256(p)!=h]
    assert not changed, changed
    return dict(status='PASS', swcs_checked=len(checked), protected_files=len(snapshot),
                exact_decision_replay=True, exact_geometry_radius_type_edges=True,
                registration_rerun=False, nn_rerun=False, vascularmd_rerun=False)


def run(input_dir, output, *, model=True):
    cache = load_cache(input_dir)
    if output.exists():
        raise FileExistsError('Derived output exists; use a new --output to preserve all prior runs: ' + str(output))
    snapshot = protection(cache, output)
    source = read_source(Path(cache.manifest['source']['raw_source']))
    for n, p in cache.points.items():
        np.testing.assert_array_equal([float(p[k]) for k in ('x','y','z','radius')], source.graph.nodes[n]['coords'])
        if int(p['swc_type']) != source.graph.nodes[n]['swc_type']:
            raise ValueError('Frozen point TYPE differs from raw source')
    output.mkdir(parents=True)
    write_json(output / 'protection_before.json', snapshot)
    summary = dict(subject='BG001', side='RMCA', rules=RULES, registration_rerun=False, nn_rerun=False,
        source=cache.provenance, components=[], exports=[], models={}, excluded_derived_components=[],
        borderline_branches=[dict(branch_id=k, component=cache.branch_component[k], p_MeVO=b['p_MeVO'], status='SEMANTIC_BORDERLINE')
            for k,b in cache.branches.items() if k in cache.branch_component and b['p_MeVO'] <= .65])
    edges = branch_edges(cache)
    lengths = {k: b['length_mm'] for k,b in cache.branches.items()}
    proximal, diameter, proximal_audit, diameter_audit = [], [], [], []
    originals, semantic, filtered = {}, {}, {}
    tasks = []
    try:
        with protect_inputs(snapshot), frozen_upstream() as disabled:
            summary['disabled_upstream_entrypoints'] = disabled
            for component in cache.components:
                old, graphs, rows, audit = refine_component(cache, component, source)
                originals[component], semantic[component], filtered[component] = old, graphs, []
                proximal.extend(rows)
                proximal_audit.extend(audit)
                for i, graph in enumerate(graphs, 1):
                    derived = f'part{component:02d}' + (f'_{i:02d}' if len(graphs)>1 else '')
                    pruned, cuts, details = prune_diameter(graph, cache.coordinates, edges, component, derived)
                    diameter.extend(cuts)
                    diameter_audit.extend(details)
                    filtered[component].append(pruned)
                    for level, g, suffix in [('semantic_refined',graph,'semantic_refined'), ('diameter_075',pruned,'d075')]:
                        if not g:
                            summary['excluded_derived_components'].append(derived)
                            continue
                        path = output / level / f'roi_{derived}_{suffix}.swc'
                        exported = export_subset(path, g, source, edges, output / 'mappings')
                        exported.update(component=component, derived_component=derived, level=level)
                        summary['exports'].append(exported)
                        tasks.append((component, level, path, g, exported))
                sem = nx.compose_all(graphs) if graphs else nx.DiGraph()
                final = nx.compose_all(filtered[component]) if filtered[component] else nx.DiGraph()
                original_keys = {edges[e] for e in sem.edges}
                final_keys = {edges[e] for e in final.edges}
                affected_keys = {edges[e] for e in set(sem.edges)-set(final.edges)}
                summary['components'].append(dict(component=component,
                    original=statistics(old,cache.coordinates,edges,lengths),
                    semantic=statistics(sem,cache.coordinates,edges,lengths),
                    diameter=statistics(final,cache.coordinates,edges,lengths),
                    removed_node_count=len(set(sem)-set(final)),
                    removed_length_mm=graph_length(sem,cache.coordinates)-graph_length(final,cache.coordinates),
                    removed_original_terminal_count=sum(sem.out_degree(n)==0 and n not in final for n in sem),
                    fully_removed_branches=sorted(original_keys-final_keys),
                    partly_removed_branches=sorted(affected_keys & final_keys),
                    sustained_cuts=sum(r['component']==component and r['reason']=='DISTAL_DIAMETER_CUT' for r in diameter)))
            summary['proximal_refinement'] = proximal
            summary['diameter_pruning'] = diameter
            table(output/'proximal_refinement.csv',proximal)
            table(output/'diameter_pruned_segments.csv',diameter)
            write_json(output/'decision_audit.json',dict(proximal=proximal_audit,diameter=diameter_audit))
            c = summary['components']
            summary['totals'] = dict(original_length_mm=sum(v['original']['length_mm'] for v in c),
                semantic_length_mm=sum(v['semantic']['length_mm'] for v in c),
                diameter_length_mm=sum(v['diameter']['length_mm'] for v in c),
                semantic_nodes=sum(v['semantic']['node_count'] for v in c),
                diameter_nodes=sum(v['diameter']['node_count'] for v in c),
                diameter_removed_nodes=sum(v['removed_node_count'] for v in c),
                diameter_removed_length_mm=sum(v['removed_length_mm'] for v in c),
                removed_original_terminals=sum(v['removed_original_terminal_count'] for v in c),
                fully_removed_branches=sum(len(v['fully_removed_branches']) for v in c),
                partly_removed_branches=sum(len(v['partly_removed_branches']) for v in c),
                sustained_cuts=sum(v['sustained_cuts'] for v in c),
                rebound_vetoes=sum(r['reason']=='DIAMETER_REBOUND' for r in diameter))
            print('ALL DERIVED SWC QC PASS', json.dumps(summary['totals']), flush=True)
            write_json(output/'refinement_data_summary.json', summary)
            export_figures(cache,originals,semantic,filtered,proximal,diameter,output/'QC')
            if model:
                # Both part03 inputs first; all SWCs have already passed QC.
                for component,level,path,graph,exported in sorted(tasks,key=lambda t:(t[0]!=3,t[0],t[1]!='semantic_refined',str(t[2]))):
                    model_dir = output/'VascularMD'/level/path.stem
                    result = model_one(path,model_dir)
                    result['single_inlet_stem'] = exported['single_inlet_stem']
                    summary['models'][path.stem] = result
                    failures = result['failures']
                    if any('branching inlet' in f['message'] for f in failures):
                        root = validate_tree(graph)
                        incoming = [k for k,b in cache.branches.items() if b['node_ids'][-1]==root]
                        if len(incoming)!=1:
                            raise ValueError('Cannot unambiguously select real upstream context')
                        context = set(cache.branches[incoming[0]]['node_ids'])-set(graph)
                        model_graph = source.graph.subgraph(set(graph)|context).copy()
                        context_path = output/'VascularMD'/'model_inputs'/(path.stem+'_with_context.swc')
                        info = export_subset(context_path,model_graph,source,edges,output/'mappings',context_nodes=context)
                        context_result = model_one(context_path,model_dir/'with_real_context')
                        context_result.update(single_inlet_stem=info['single_inlet_stem'],context_branch=incoming[0],
                                              context_export=info,not_part_of_semantic_roi=True)
                        summary['models'][context_path.stem] = context_result
                    write_json(model_dir/'model_export_summary.json',result)
            write_json(output/'refinement_summary.json',summary)
            (output/'refinement_report.md').write_text(report_text(summary),encoding='utf-8')
    finally:
        changed = [p for p,h in snapshot.items() if not Path(p).is_file() or sha256(p)!=h]
        write_json(output/'protection_verification.json',dict(protected_files=len(snapshot),changed=changed,
                   all_unchanged=not changed,registration_rerun=False,nn_rerun=False))
        if changed:
            raise ValueError('Protected source changed: '+str(changed))
    print('COMPLETE', output, flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--skip-model',action='store_true',help='Explicit data-only run, still exports and validates both SWC layers')
    parser.add_argument('--verify',action='store_true',help='Read-only replay of cuts and verification of existing outputs; no modeling')
    args = parser.parse_args()
    output = (args.output or args.input/'refined_roi').resolve()
    if args.verify:
        print(json.dumps(verify_results(args.input.resolve(),output),indent=2))
    else:
        run(args.input.resolve(), output, model=not args.skip_model)


if __name__ == '__main__':
    main()
