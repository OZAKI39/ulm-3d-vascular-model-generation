#!/usr/bin/env python3
"""Export the frozen BG001 RMCA boundary evidence without recomputing production."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
import shutil
import tempfile
import zipfile

import networkx as nx

from vascular_processing.boundary_review_export import (load_cache,build_evidence,separation_pairs,
    component_summary,local_keys,prohibit_recomputation,protect_inputs,ReviewSourceMismatch,TYPES,write_csv)
from vascular_processing.topbrain_qc import sha256,write_json


DEFINITIONS={
 'PROXIMAL_M1_TO_MEVO':'上游 M1 分支与下游 MeVO 分支之间的候选入口。',
 'DISTAL_MEVO_TO_UNKNOWN':'MeVO 结束并进入 UNKNOWN 的候选远端。',
 'DISTAL_MEVO_TO_TERMINAL':'MeVO 分支在原始 SWC 终端结束；不等同于已确认的 M3 末端。',
 'UNKNOWN_TO_MEVO':'UNKNOWN 下游进入 MeVO；需审查支持空缺及入口依据。',
 'MEVO_TO_M1_REVERSE':'MeVO 下游重新出现 M1 的逆序候选；最高优先级。',
 'MEVO_COMPONENT_ROOT':'每个严格连通组件的首层 MeVO 分支；分叉根可对应多条首层分支。',
 'DETACHED_MEVO_COMPONENT':'组件入口没有直接上游 M1；不表示原始血管断裂，也不判定语义错误。'}


def copy_evidence(cache,out):
    evidence=out/'evidence';evidence.mkdir()
    for name in ['roi_manifest.json','point_predictions.csv','branch_predictions.csv','labeled_tree.vtp',
        'node_mapping.csv','branch_mapping.csv','donor_registrations.csv','donor_registrations.json','point_support.npz']:
        p=cache.directory/name
        if p.exists():shutil.copy2(p,evidence/name)
    shutil.copy2(cache.directory.parent/'input_provenance.json',evidence/'input_provenance.json')
    for key in ['raw_source','color_source']:
        path=Path(cache.manifest['source'][key]);shutil.copy2(path,evidence/path.name)
    for p in (cache.directory/'ROI').iterdir():
        if p.is_file():
            (evidence/'ROI').mkdir(exist_ok=True);shutil.copy2(p,evidence/'ROI'/p.name)
    for component in cache.components.values():
        source=Path(component['VascularMD_report']);folder=evidence/'VascularMD_QC'/f'roi_part{component["component"]:02d}';folder.mkdir(parents=True)
        shutil.copy2(source,folder/source.name)
        for log in source.parent.glob('*_run.log'):shutil.copy2(log,folder/log.name)
    write_json(evidence/'source_hash_verification.json',dict(status='PASS',provenance=cache.provenance,
        source_files=[dict(file=p,sha256=h) for p,h in cache.snapshot.items()],
        omitted_large_artifacts=['TopBrain NIfTI','full donor semantic clouds','native VascularMD surfaces/STL'],
        note='All necessary review geometry is included as original SWC and compact line VTP. Original large-surface paths are provenance only.'))


def branch_table(cache,events):
    boundary_related={k for e in events for k in local_keys(cache,e)}
    return [dict(**{key:b[key] for key in ['branch_id','parent_branch','label','p_M1','p_MeVO','known_fraction','confidence',
        'length_mm','radius_median','valid_donor_count','support_distance']},roi_component=cache.branch_component.get(k),
        is_boundary_related=k in boundary_related) for k,b in cache.branches.items()]


def write_topology(cache,events,out):
    identifiers={k:[] for k in cache.branches}
    for e in events:
        key=e['downstream_branch_id'] or e['upstream_branch_id']
        if key in identifiers:identifiers[key].append(e['boundary_id'])
    lines=['BG001 RMCA root | native branch-parent direction | world units: mm']
    def walk(key,prefix,last):
        b=cache.branches[key];p='null' if b['p_MeVO'] is None else f'{b["p_MeVO"]:.4f}'
        lines.append(prefix+('└── ' if last else '├── ')+f'branch {key} [{b["label"]}] pMeVO={p}'+
            (' <BOUNDARY: '+', '.join(identifiers[key])+'>' if identifiers[key] else ''))
        children=sorted(cache.branch_graph.successors(key))
        for i,child in enumerate(children):walk(child,prefix+('    ' if last else '│   '),i==len(children)-1)
    roots=[k for k in cache.branch_graph if cache.branch_graph.in_degree(k)==0]
    for i,key in enumerate(roots):walk(key,'',i==len(roots)-1)
    (out/'tables/BG001_RMCA_topology.txt').write_text('\n'.join(lines)+'\n')
    write_json(out/'tables/BG001_RMCA_topology.json',dict(root_original_swc_id=cache.root,root_branches=roots,
        nodes=[dict(b,roi_component=cache.branch_component.get(k),boundary_ids=identifiers[k]) for k,b in cache.branches.items()],
        edges=[dict(parent=a,child=b) for a,b in cache.branch_graph.edges]))
    write_json(out/'tables/rmca_root.json',dict(root_node=cache.root,root_branch=roots,coordinate=cache.coordinates[cache.root][:3],
        radius=float(cache.coordinates[cache.root][3]),units='mm',major_tree_TYPE=4,
        original_standardized_swc_TYPE=int(cache.points[cache.root]['swc_type']),
        type_note='Major-tree TYPE 4 is companion ColorCoded metadata; the standardized source TYPE is preserved',mapping_provenance=cache.provenance))


def markdown_docs(cache,events,components,pairs,summary,out):
    definitions='\n'.join(f'- **{key}**：{value}' for key,value in DEFINITIONS.items())
    readme=f'''# BG001 RMCA MeVO 候选边界：外部审查入口

本包包含一个锁定版本的 BraVa BG001 右侧大脑中动脉及其 TopBrain 引导的语义候选。MeVO = M2 OR M3；M1 为近端类别，MeVO 为合并后的 M2/M3 候选，UNKNOWN 表示当前支持不足或分支判别未满足固定阈值。标签不是人工 ground truth。本次只整理已经冻结的结果，没有重新配准、投票、聚合分支、提取 ROI 或建模。

TYPE 4 = RMCA 采用 **BG001-specific DERIVED_SOURCE_SUPPORTED** 约定：由官方 BG001 具名长度表推断，并由用户明确选定。仅绑定本包原始和 ColorCoded SWC 的哈希，不推广至其他 BraVa。生产记录中的历史左右命名冲突仍完整保留在 evidence；本次审查状态命名不表示获得了新的官方直接 TYPE 字典。

## 建议阅读顺序

1. 阅读 `EXTERNAL_REVIEW_SUMMARY.md`，再查看 `global/` 中的三类总览。边界总览另提供相差 90° 的第二视图。
2. 优先阅读 `HIGH_PRIORITY_REVIEW.md`。每项链接到固定四面板图片和机器可读证据。
3. 图 A：完整 RMCA 为淡灰色，当前组件突出，红色圆形 marker 为当前边界。图 B：最多两条上游分支及两层下游分支；图 C：同一区域第二视角。图 D：直接来自缓存的概率、支持距离、半径、长度与供体证据。
4. PROXIMAL 表示靠近原生 RMCA 根的上游侧，DISTAL 表示沿原生 parent→child 的下游侧。箭头遵从 SWC 父子方向，不能等同于已经测量的生理血流方向。
5. 对照 `tables/` 的 CSV 和缩进拓扑。必要时使用单边界的 metrics/local_topology JSON；VTP 仅是可选的第二级证据。

每张图有固定颜色图例：M1 橙色，MeVO 青色，UNKNOWN 紫色，当前边界红色，背景 context 灰色。组件总览另以五种颜色区分 part01–part05。中心线按原样绘制，线宽是显示用宽度；真实半径以 CSV/JSON 数值为准，不以 PNG 线宽测量。

## 边界类型

{definitions}

一个物理分叉点可同时具有入口、组件根等多个事件；分叉根的每条首层 MeVO 分支也单独保留事件。ID 由 subject、side、上游分支、下游分支和类型确定，与 CSV 行号无关。因此事件数不等于空间位置数。总览 marker 后的 `+` 表示同一位置有多个事件，旁侧索引列出同位置的完整短 ID；每个事件都有独立图片。

## 数值和状态的边界

概率、known fraction、confidence、分支支持距离和有效 donor 数全部来自生产缓存。组件平均概率仅是对这些既有分支指标按已存长度作描述性汇总，不改变任何分支标签。radius ratio = 下游分支原始样本半径中位数 / 上游对应中位数。根距离及最近分叉距离均沿原生中心线弧长测量；最近分叉允许包含边界点自身。

供体 before/after 证据读取共享节点前、后最近的**不同原始 SWC 样本**，同时给出节点编号；不是重新发起 NN 查询。它可能不同于全分支的平均支持。终端不存在下游分支，相关字段为 null，CSV 对应空单元；这不是零概率。任何缺字段会明确列为 FIELD_NOT_AVAILABLE，不猜值。

所有事件的 review_status 保持 **UNREVIEWED**。高优先级仅表示应优先审查，不表示 ACCEPT、REJECT 或已确定的 M2 起点。

**VascularMD entry failure does NOT imply semantic boundary failure.** part03 是有效保留的语义候选区域，其建模入口检查失败与语义边界是否合理完全分开。本包不会把建模失败自动转为语义否定。所有几何仍为 ANATOMICAL_TRANSFER_CANDIDATE，制造状态为 MANUFACTURING_NOT_VALIDATED。

## 文件结构与离线使用

```text
README_FOR_EXTERNAL_REVIEW.md / EXTERNAL_REVIEW_SUMMARY.md / HIGH_PRIORITY_REVIEW.md
review_manifest.json / bundle_generation_summary.json / bundle_generation.log
tables/       边界、组件、分支 CSV；55 分支拓扑 TXT/JSON；根与组件分离原因
global/       标签、组件、边界总览及边界第二视角
components/   五个组件总览
boundaries/<boundary_id>/  review PNG、metrics JSON、local_topology JSON、VTP+显示说明
evidence/     两个 BG001 源 SWC、原始预测和映射、小型支持 NPZ、严格/modelable SWC、原生 QC 与日志
```

生产表中的绝对路径用于溯源；查看本包无需这些路径存在。严格/modelable SWC 的 review 表使用包内相对路径。大型 VascularMD surface/STL 仅记录原保存路径，未复制；其状态由包内 QC 和日志说明。未包含 TopBrain NIfTI、供体点云、其他 BraVa 个体或完整 VascularMD 缓存。PNG、CSV、JSON 和源中心线足够支持第一轮边界审查。

`review_manifest.json` 对每个普通文件记录 SHA256、用途、优先级和关联事件/组件。manifest 本身使用“将自己的 SHA256 字段设为 null 后的 canonical JSON”的摘要规则，避免自引用哈希不可能固定的问题；完整 ZIP 摘要位于包旁的 `.sha256` 文件。ZIP 内路径均为相对路径。
'''
    (out/'README_FOR_EXTERNAL_REVIEW.md').write_text(readme,encoding='utf-8')
    counts=summary['boundary_counts_by_type']
    count_table='\n'.join(f'| {kind} | {counts[kind]} |' for kind in TYPES)
    component_table='\n'.join(f'| {r["component_id"]} | {r["node_count"]} / {r["branch_count"]} | {r["proximal_branch_id"]} | {r["proximal_parent_label"]} | {r["separation_reason"]} |' for r in components)
    boundary_table='\n'.join(f'| [{e["boundary_id"].split("_")[-1]}](boundaries/{e["boundary_id"]}/{e["boundary_id"]}_review.png) | {e["boundary_type"]} | {e["roi_component"]:02d} | {e["upstream_branch_id"] or "—"} → {e["downstream_branch_id"] or "TERMINAL"} |' for e in events)
    text=f'''# BG001 RMCA 候选边界审查摘要

## A. 来源与范围

审查对象为已冻结的 BG001 RMCA、TYPE 4。映射状态 DERIVED_SOURCE_SUPPORTED，仅适用于当前指定文件。raw SHA256：`{RAW(cache)}`；ColorCoded SHA256：`{cache.manifest['source']['color_sha256']}`。分支 CSV SHA256：`{cache.provenance['branch_predictions_sha256']}`。完整 ROI manifest 源绑定及下载版本说明见包内 evidence。

本次导出只读取生产结果。registration、ranking、NN、UNKNOWN 校准、分支聚合、ROI 提取、VascularMD 建模及表面生成均未执行；真实导出时相应入口被禁止调用断言保护。所有审查状态仍为 UNREVIEWED。

## B. RMCA 总体统计

共 797 个原始节点、55 条原生拓扑分支：M1 8 条，MeVO 44 条，UNKNOWN 3 条。严格 MeVO 为 5 个组件，合计 687 个原始节点、44 条分支。完整分支保留共享端点，因此与点级 MeVO 的 678 个节点数量不同。供体为 MRA024、013、015、022、011；既有支持阈值为 agreement≥0.60、有效 donor≥3、归一化距离≤12.7949772843。

## C. 候选事件计数

共 **{len(events)} 个事件**，对应 **{summary['boundary_location_count']} 个原始 SWC 位置**。同一组件分叉根可有多条首层分支，也可同时产生不同类型事件。

| 类型 | 数量 |
|---|---:|
{count_table}

## D. 五个组件与分离原因

| 组件 | 节点/分支 | 首层 MeVO 分支 | 直接上游标签 | 与其他组件的分离原因集合 |
|---|---:|---|---|---|
{component_table}

分离原因根据既有 branch-parent 路径描述，未修改任何连接。完整两两关系、连接路径和各路径标签见 `tables/component_separation.json`。UNKNOWN_GAP 指相关组件之间的原生路径含 UNKNOWN；DIFFERENT_PROXIMAL_PARENT 指它们具有不同的直接近端父分支；其余原因的定义见该 JSON。

## E. 全部边界简表

短 ID 链接到独立四面板图；完整 ID、原始节点号和全部数值在边界 CSV 中。

| 短 ID | 类型 | part | 上游 → 下游 |
|---|---|---:|---|
{boundary_table}

## F. 高优先级审查

共有 **{summary['high_priority_count']} 项 HIGH_PRIORITY_REVIEW**。包括全部近端 M1→MeVO、UNKNOWN→MeVO、每个组件首层分支，以及存在时的逆序。支持较低的 MeVO/UNKNOWN 转换按缓存 known fraction、confidence 排在同类前部；不因排序改变其标签或审查状态。对应图片索引见 `HIGH_PRIORITY_REVIEW.md`。

特别应检查 part03 的两个 UNKNOWN→MeVO 入口及其共同分叉位置。此处只提出审查问题：供体支持、相邻节点投票与组件入口的关系是否足够清楚？没有自动认定该候选正确或错误。

## G. VascularMD 状态

part01、02、04、05 已在此前完成原生建模及 smooth SWC/VTK/STL。part03 的原生入口检查失败：其根为分叉且没有符合既有规则的直接上游 M1 context。**VascularMD entry failure does NOT imply semantic boundary failure.** 本次不重跑模型，不改变 context，不产生新血管表面。QC 和原运行日志已独立复制入包。

## H. 已知限制

这些是 TopBrain-informed candidate labels，不是 BraVa 人工分段真值；Pilot 指标不能当作 BG001 准确率。既有 OOD 支持阈值未被本次重新优化。类型4的左右命名是 BG001 特定来源推断。PNG 是中心线证据图，不能用于精确量取管腔外表面；半径以原始数值为准。所有方向箭头只代表 native parent→child。终端也不意味着已确认的 M3 末端。

模型与制造适用性均不代替语义审查。所有真实事件均为 UNREVIEWED。缺数据或渲染错误的完整列表位于 `bundle_generation_summary.json`；本次 warnings：{'; '.join(summary['warnings']) or '无'}。Errors：{'; '.join(summary['errors']) or '无'}。
'''
    (out/'EXTERNAL_REVIEW_SUMMARY.md').write_text(text,encoding='utf-8')
    high=[e for e in events if e['priority']!='MEDIUM']
    high.sort(key=lambda e:({'MEVO_TO_M1_REVERSE':0,'PROXIMAL_M1_TO_MEVO':1,'UNKNOWN_TO_MEVO':2,'DISTAL_MEVO_TO_UNKNOWN':2,'MEVO_COMPONENT_ROOT':3}.get(e['boundary_type'],4),
        min([e[s+'_known_fraction'] for s in ['upstream','downstream'] if e[s+'_known_fraction'] is not None] or [1]),
        min([e[s+'_confidence'] for s in ['upstream','downstream'] if e[s+'_confidence'] is not None] or [1]),
        -max([e[s+'_median_support_distance'] for s in ['upstream','downstream'] if e[s+'_median_support_distance'] is not None] or [0]),e['boundary_id']))
    (out/'HIGH_PRIORITY_REVIEW.md').write_text('# HIGH_PRIORITY_REVIEW\n\n仅为审查顺序，全部 UNREVIEWED。\n\n'+
        '\n'.join(f'- [{e["boundary_id"]}](boundaries/{e["boundary_id"]}/{e["boundary_id"]}_review.png)：{e["boundary_type"]}，part{e["roi_component"]:02d}，{e["priority"]}。' for e in high)+'\n',encoding='utf-8')


def RAW(cache):return cache.manifest['source']['raw_sha256']


def make_manifest(out,events):
    by_id={e['boundary_id']:e for e in events};files=[]
    for p in sorted(out.rglob('*')):
        if not p.is_file() or p.name=='review_manifest.json':continue
        relative=str(p.relative_to(out));boundary=next((part for part in p.parts if part in by_id),None)
        component=by_id[boundary]['roi_component'] if boundary else next((part[:10] for part in p.parts if part.startswith('roi_part')),None)
        files.append(dict(file=relative,type=p.suffix.lstrip('.') or 'text',boundary_id=boundary,component=component,
            SHA256=sha256(p),description=('Boundary evidence: ' if boundary else 'Review package: ')+p.name,
            priority=by_id[boundary]['priority'] if boundary else 'PACKAGE',size_bytes=p.stat().st_size))
    files.append(dict(file='review_manifest.json',type='json',boundary_id=None,component=None,SHA256=None,
        description='Package inventory; SHA256 uses canonical JSON with this self SHA256 field set to null',priority='PACKAGE'))
    manifest=dict(schema='bg001-rmca-external-review-v1',subject='BG001',side='RMCA',
        self_hash_rule='SHA256 of UTF-8 json.dumps(manifest, sort_keys=True, separators=(comma,colon), ensure_ascii=False), with the self-entry SHA256 set to null',files=files)
    canonical=json.dumps(manifest,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()
    files[-1]['SHA256']=hashlib.sha256(canonical).hexdigest()
    write_json(out/'review_manifest.json',manifest)


def pack(out,path):
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for file in sorted(out.rglob('*')):
            if file.is_file():archive.write(file,Path('BG001_RMCA_external_review')/file.relative_to(out))


def generate(input_dir,output_dir):
    from vascular_processing.boundary_review_figures import boundary_figure,overview,write_local_vtp
    # No output directory or ZIP is created until all mandatory hashes match.
    with prohibit_recomputation() as guards:
        cache=load_cache(input_dir)
        output=Path(output_dir).resolve();archive=output.parent/'BG001_RMCA_external_review_bundle.zip'
        if output==cache.directory or output in [Path(p) for p in cache.snapshot]:raise ValueError('Output must be a separate review directory')
        if output.exists() or archive.exists():
            # Re-running the documented command verifies and reuses this exact
            # review package. Never overwrite an edited or incomplete review.
            evidence=json.loads((output/'evidence/source_hash_verification.json').read_text())
            if evidence['provenance']!=cache.provenance:raise ReviewSourceMismatch('Existing review belongs to different production evidence')
            manifest=json.loads((output/'review_manifest.json').read_text())
            for item in manifest['files']:
                if item['file']!='review_manifest.json' and sha256(output/item['file'])!=item['SHA256']:
                    raise ValueError('REVIEW_BUNDLE_INTEGRITY_FAILED: existing review was edited; preserve it and choose another output directory')
            own=next(item for item in manifest['files'] if item['file']=='review_manifest.json')
            expected=own['SHA256'];own['SHA256']=None
            if hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()!=expected:
                raise ValueError('REVIEW_BUNDLE_INTEGRITY_FAILED: manifest self hash')
            if sha256(archive)!=archive.with_suffix('.zip.sha256').read_text().split()[0]:raise ValueError('REVIEW_BUNDLE_INTEGRITY_FAILED: ZIP hash')
            cache.verify_unchanged()
            print('Existing review verified; no figures or production stages recomputed.\nReview bundle ready:\n'+str(archive),flush=True)
            return json.loads((output/'bundle_generation_summary.json').read_text())
        output.parent.mkdir(parents=True,exist_ok=True)
        with protect_inputs(cache.snapshot),tempfile.TemporaryDirectory(prefix='.review_build_',dir=output.parent) as temporary:
            out=Path(temporary)/'manual_review';out.mkdir();messages=[]
            def log(message):
                line=datetime.now(timezone.utc).isoformat(timespec='seconds')+' '+message
                messages.append(line);print(line,flush=True)
            log('INPUT '+str(cache.directory));log(f'HASH_VERIFICATION PASS: {len(cache.snapshot)} immutable inputs')
            for path,digest in cache.snapshot.items():messages.append(f'INPUT_SHA256 {digest} {path}')
            events,details,local=build_evidence(cache);pairs=separation_pairs(cache);components=component_summary(cache,events,pairs)
            for folder in ['tables','global','components','boundaries']:(out/folder).mkdir()
            copy_evidence(cache,out)
            write_csv(out/'tables/BG001_RMCA_boundary_candidates.csv',events)
            write_csv(out/'tables/BG001_RMCA_component_summary.csv',components)
            write_csv(out/'tables/BG001_RMCA_branch_predictions_review.csv',branch_table(cache,events))
            write_json(out/'tables/component_separation.json',dict(pairs=pairs,policy='All ten component pairs; shortest existing branch-parent path; no reconnection',
                allowed_reasons=['UNKNOWN_GAP','DIFFERENT_PROXIMAL_PARENT','TOPOLOGICALLY_SEPARATE_DAUGHTERS','NO_DIRECT_CONNECTION','OTHER']))
            write_topology(cache,events,out)
            errors=[];warnings=list(cache.warnings)
            log(f'CANDIDATE_DETECTION {len(events)} events, {len({e["boundary_node_id"] for e in events})} physical nodes')
            for i,e in enumerate(events,1):
                key=e['boundary_id'];folder=out/'boundaries'/key;folder.mkdir()
                write_json(folder/'metrics.json',details[key]);write_json(folder/'local_topology.json',local[key])
                write_local_vtp(cache,e,local[key],folder/(key+'_review.vtp'))
                try:boundary_figure(cache,e,details[key],local[key],folder/(key+'_review.png'))
                except Exception as exc:
                    errors.append(f'BOUNDARY_RENDER_FAILED: {key}: {type(exc).__name__}: {exc}');log(errors[-1])
                for side in details[key]['donor_vote_breakdown'].values():
                    if side.get('warning'):warnings.append(f'{side["warning"]}: {key}, node {side["node_id"]}')
                log(f'BOUNDARY {i}/{len(events)} {key} '+('PNG_OK' if (folder/(key+'_review.png')).is_file() else 'PNG_FAILED'))
            for prefix,kind in [('00','labels'),('01','components'),('02','boundaries')]:
                overview(cache,events,out/'global'/f'{prefix}_BG001_RMCA_{kind}.png',kind)
            overview(cache,events,out/'global/02_BG001_RMCA_boundaries_viewB.png','boundaries',azim=25)
            for component in cache.components:overview(cache,events,out/'components'/f'roi_part{component:02d}_overview.png','component',component=component)
            cache.verify_unchanged();log('READ_ONLY PASS; registration/NN/calibration/aggregation/extraction/VascularMD guards NOT CALLED')
            label_counts=Counter(b['label'] for b in cache.branches.values())
            summary=dict(subject='BG001',side='RMCA',source_hashes=cache.manifest['source'],mapping_status='DERIVED_SOURCE_SUPPORTED',
                branch_count=len(cache.branches),M1_branch_count=label_counts['M1'],MeVO_branch_count=label_counts['MeVO'],UNKNOWN_branch_count=label_counts['UNKNOWN'],
                component_count=len(components),boundary_count=len(events),boundary_location_count=len({e['boundary_node_id'] for e in events}),
                boundary_counts_by_type={kind:sum(e['boundary_type']==kind for e in events) for kind in TYPES},
                high_priority_count=sum(e['priority']!='MEDIUM' for e in events),
                PNG_count=len(list(out.rglob('*.png'))),VTM_count=0,VTP_count=len(list(out.rglob('*.vtp'))),
                boundary_VTP_count=len(list((out/'boundaries').rglob('*.vtp'))),cached_context_VTP_count=1,
                total_bundle_size=0,size_units='bytes',size_measurement='Exact final ZIP size',
                warnings=sorted(set(warnings)),errors=errors,review_status='UNREVIEWED',
                all_boundaries_unreviewed=all(e['review_status']=='UNREVIEWED' for e in events),
                recomputation_guarded_entry_points=guards,recomputation_count=0,protected_input_count=len(cache.snapshot),
                read_only_verified=True,s1_2_sha256=sha256(ROOT/'s1-2_swc_roi_generate_human.py'),components=components)
            markdown_docs(cache,events,components,pairs,summary,out)
            summary['file_count']=len([p for p in out.rglob('*') if p.is_file()])+3
            for warning in summary['warnings']:log('WARNING '+warning)
            # Uncompressed members with reserved JSON whitespace let the exact
            # integer ZIP size be filled in without a self-referential size cycle.
            write_json(out/'bundle_generation_summary.json',summary)
            summary_path=out/'bundle_generation_summary.json'
            summary_size=summary_path.stat().st_size+64
            summary_path.write_bytes(summary_path.read_bytes().ljust(summary_size,b' '))
            log('PACKAGING: all manual_review files; no NIfTI, donor clouds, or native surfaces')
            messages.append('FINAL_ZIP_SIZE_BYTES=00000000000000000000')
            (out/'bundle_generation.log').write_text('\n'.join(messages)+'\n')
            make_manifest(out,events)
            temp_zip=Path(temporary)/archive.name
            def pack_fixed():
                with zipfile.ZipFile(temp_zip,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
                    for p in sorted(out.rglob('*')):
                        if p.is_file():
                            compression=zipfile.ZIP_STORED if p.name in {'bundle_generation_summary.json','bundle_generation.log','review_manifest.json'} else zipfile.ZIP_DEFLATED
                            z.write(p,Path('BG001_RMCA_external_review')/p.relative_to(out),compress_type=compression)
            pack_fixed();size=temp_zip.stat().st_size
            summary['total_bundle_size']=size;write_json(summary_path,summary)
            summary_path.write_bytes(summary_path.read_bytes().ljust(summary_size,b' '))
            messages[-1]=f'FINAL_ZIP_SIZE_BYTES={size:020d}'
            (out/'bundle_generation.log').write_text('\n'.join(messages)+'\n');make_manifest(out,events);pack_fixed()
            assert temp_zip.stat().st_size==size,'Unexpected package size feedback'
            if size>=100*1024**2:raise ValueError('Review ZIP exceeds 100 MiB; review static image compression')
            with zipfile.ZipFile(temp_zip) as z:
                if z.testzip() is not None:raise ValueError('ZIP CRC verification failed')
            cache.verify_unchanged()
            shutil.move(str(out),output);shutil.move(str(temp_zip),archive)
        digest=sha256(archive);archive.with_suffix('.zip.sha256').write_text(digest+'  '+archive.name+'\n')
        print('Review bundle ready:\n'+str(archive),flush=True)
        return summary


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA')
    p.add_argument('--output-dir',type=Path,default=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/manual_review')
    args=p.parse_args(argv)
    try:
        summary=generate(args.input_dir,args.output_dir)
        return 0 if not summary['errors'] else 2
    except Exception as exc:
        print(str(exc),file=sys.stderr,flush=True);return 2


if __name__=='__main__':raise SystemExit(main())
