"""Instrument the frozen V1 source; distinguish replay from historical logs."""
import json
from collections import Counter
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save,table,digest,transform
from . import support_custom_graph as original
from .support_custom_pipeline import bottom_mesh
from .support_custom_qc import self_printability

def run(out,inputs):
    legacy=out.parent/'legacy_explicit_support_v1';report=legacy/'vascular_support_report.md';text=report.read_text();summary=json.loads((legacy/'vascular_support_summary.json').read_text())
    path=legacy/'source_snapshot/vascular_processing/support_custom_graph.py';source=path.read_text()
    source=source.replace('nodes=[];edges=[];pieces=[];tips=[];roots=[];feet=[];rejected=[];repairs=[];active=[]','nodes=[];edges=[];pieces=[];tips=[];roots=[];feet=[];rejected=[];repairs=[];active=[];merge_audit=[]')
    old="if m[2]<floor+4 or np.linalg.norm(a-b)>cfg.get('merge_max_distance_mm',8.):continue"
    new="""if m[2]<floor+4:
                merge_audit.append(dict(node_a=i,node_b=j,reason='root_restriction',merge_xyz=m.tolist()));continue
            if np.linalg.norm(a-b)>cfg.get('merge_max_distance_mm',8.):
                merge_audit.append(dict(node_a=i,node_b=j,reason='neighbourhood_cutoff',merge_xyz=m.tolist()));continue"""
    if old not in source:raise ValueError('V1_SNAPSHOT_CHANGED_CANNOT_INSTRUMENT')
    source=source.replace(old,new)
    target='if all(legal(v) for v in arms):choices.append((float(np.linalg.norm(a-m)),j,m,d,da,db))'
    replacement="""valid=all(legal(v) for v in arms)
            merge_audit.append(dict(node_a=i,node_b=j,reason='legal_pair' if valid else 'collision_or_xy_gap',merge_xyz=m.tolist()))
            if valid:choices.append((float(np.linalg.norm(a-m)),j,m,d,da,db))"""
    if target not in source:raise ValueError('V1_MERGE_PROBE_SOURCE_MISMATCH')
    source=source.replace(target,replacement).replace("graph=dict(nodes=nodes,edges=retained_edges","graph=dict(merge_audit=merge_audit,nodes=nodes,edges=retained_edges")
    namespace=dict(original.__dict__);exec(compile(source,str(path), 'exec'),namespace)
    cfg=json.loads((legacy/'resolved_config.json').read_text());pr=summary['profile'];contacts=json.loads((legacy/'contact_plan.json').read_text())['contacts']
    v=trimesh.load_mesh(inputs['core']);v.apply_transform(inputs['print_transform']);box=trimesh.load_mesh(inputs['box']);box.apply_transform(inputs['print_transform']);inner=transform(inputs['inner_bounds'],inputs['print_transform'])
    selected=summary['selected'];support,g,feet=namespace['build'](contacts,v,box,inner,pr,cfg,selected['tip_width_mm'],selected['neck_width_mm'],selected['pillar_diameter_mm'])
    check=self_printability(support,bottom_mesh(box,inner),pr);counts=Counter(r['reason'] for r in g['merge_audit']);table(out/'v1_merge_replay.csv',g['merge_audit']);save(out/'v1_merge_replay.json',dict(source_sha256=digest(path),historical_report_sha256=digest(report),source_report_excerpt=text[:600],replay_not_original_runtime_log=True,counts=dict(counts),self_printability=check,initial_merge_nodes=len(g['merge_nodes']),final_historical_merge_nodes=0))
    (out/'support_v1_postmortem.md').write_text(f'''# V1 支撑拓扑复盘

本报告自动读取已冻结的 V1 报告和摘要。V1 检测到 {summary['detection']['critical_intervals']} 个关键区间，最终有 {selected['planned_tip_count']} 个尖端、40 个根部、0 个合并节点，危险采样覆盖 {selected['planned_coverage']['percent']:.2f}%。

历史程序没有逐次保存全部合并尝试。本次对冻结的 V1 源码做了仅增加记录的确定性重放，**以下统计是重放结果，不冒充原运行日志**。源码和原报告 SHA256、每一个节点对的合并位置及原因见 v1_merge_replay.json / v1_merge_replay.csv。

|重放筛选结果|次数|
|---|---:|
|根部高度限制|{counts['root_restriction']}|
|距离超出 V1 合并邻域|{counts['neighbourhood_cutoff']}|
|碰撞或 XY 间隙不足|{counts['collision_or_xy_gap']}|
|几何允许的配对|{counts['legal_pair']}|

初始图具有 {len(g['merge_nodes'])} 个合并节点；其自身逐层检查出现 {check['unsupported_regions']} 个欠支撑区域。V1 随后的修复会把整张图重新构建为独立立柱，最终合并数因此降至 0。它没有逐个修复合并节点的高度、横向位置或桥接结构，这是全局退化的直接代码原因。

分支角度只参与原几何稳定性标记，并非上述合并筛选的独立拒绝门槛；不能编造“因角度拒绝”的历史次数。短桥不在 V1 正式拓扑回退路线中，属于未尝试，不是切片验证失败。移除发生在最后检查，未参与合并选择；V1 最终仍有 {selected['removal']['trapped_support_count']} 个被困组件。

V2 保留检测器，允许改变打印重力方向、使用盒外打印板落根，采用顶垫和量化矩形构件，并为合并逐级修复、分段颈及拆除顺序保留逐项记录。V1 仅作为 legacy_explicit_support_v1 对照，不进入新生产选择。
''')
    return dict(counts)
