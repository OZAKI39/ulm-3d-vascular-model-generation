#!/usr/bin/env python3
import collections, hashlib, json, os, re, shutil, subprocess
from pathlib import Path
from export_flow_projects import AUDIT, BASE, CE, sha

patterns={
 'private_key':re.compile(rb'-----BEGIN (?:OPENSSH|RSA|EC|DSA|PGP|ENCRYPTED)? ?PRIVATE KEY-----'),
 'github_token':re.compile(rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{35,})\b'),
 'aws_key':re.compile(rb'\bAKIA[0-9A-Z]{16}\b'),
 'openai_key':re.compile(rb'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{35,}\b'),
 'literal_credential':re.compile(rb'(?im)["\x27]?(?:api_key|access_token|auth_token|password|client_secret)["\x27]?\s*[:=]\s*["\x27]([A-Za-z0-9_+/=-]{24,})["\x27]')
}

for project in ['brava','mouse']:
 root=BASE/'temp_storage'/f'{project}_flow_github_sync_20261022';meta=root/'sync_20261022';meta.mkdir(exist_ok=True)
 inc=json.loads((AUDIT/f'{project.upper()}_INCLUDED.json').read_text());exc=json.loads((AUDIT/f'{project.upper()}_EXCLUDED.json').read_text());links=json.loads((AUDIT/f'{project.upper()}_LINKS.json').read_text())
 adaptations=[]
 if project=='brava':
  for source,destination in [(AUDIT/'remote_environment','deployment/current_server'),(BASE/'server_migration/brava_new_server_20261009/fixture','deployment/newserver_20261009/fixture')]:
   for p in source.rglob('*'):
    if not p.is_file() or p.is_symlink():continue
    dst=destination+'/'+str(p.relative_to(source));out=root/dst;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
    inc.append({'path':dst,'source':str(p),'size':p.stat().st_size,'sha256':sha(p),'mode':'100755' if p.stat().st_mode&0o111 else '100644'})
  for name in ['REMOTE_BRAVA_INVENTORY.json','REMOTE_LOCAL_COMPARISON.json']:
   shutil.copy2(AUDIT/name,meta/name)
 for row in inc:
  p=root/row['path']
  if p.name=='.gitattributes' and 'filter=lfs' in p.read_text():
   original=meta/'upstream_attributes'/row['path'];original.parent.mkdir(parents=True,exist_ok=True)
   shutil.copy2(p,original.with_name('gitattributes.original.txt'))
   p.write_text('# Snapshot packaging only; upstream rules recorded in sync_20261022/upstream_attributes.\n* -text -filter !diff !merge\n')
   adaptations.append({'path':row['path'],'source_sha256':row['sha256'],'published_sha256':sha(p),'reason':'disable inherited LFS for ordinary Git snapshot; scientific bytes unchanged'})
   row['source_sha256']=row['sha256'];row['sha256']=sha(p);row['size']=p.stat().st_size
 (root/'.gitattributes').write_text('# Preserve exact scientific bytes; this snapshot uses ordinary Git.\n* -text -filter !diff !merge\n')
 (root/'.gitignore').write_text('__pycache__/\n.pytest_cache/\n.venv/\n*.pyc\n')
 statuses=[]
 for case in json.loads((CE/'CASE_INDEX.json').read_text()):
  # Index structure is inspected below if it is not the expected list.
  if isinstance(case,dict):
   name=case.get('name','')
   if (project=='brava' and name.startswith('brava')) or (project=='mouse' and (name.startswith('mean-') or name.startswith('validation_cases/') or name=='fem_reference_stage_q')):statuses.append(case)
 (meta/'CASE_INDEX.json').write_text(json.dumps(statuses,indent=2,ensure_ascii=False)+'\n')
 for name,value in [('INCLUDED',inc),('EXCLUDED',exc),('LINKS',links),('PACKAGING_ADAPTATIONS',adaptations)]:
  (meta/f'{name}.json').write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
 for name in ['export_flow_projects.py','finalize_snapshots.py']:shutil.copy2(AUDIT/name,meta/name)
 summary={'project':project,'source_files':len(inc),'source_bytes':sum(x['size'] for x in inc),'symlinks':len(links),'excluded_entries':len(exc),'excluded_bytes':sum(x.get('size',0) for x in exc),'new_scientific_computations':0,'original_sources_modified':False,'lfs':False}
 (meta/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 branch=f'sync/{"brava-flow" if project=="brava" else "mouse-vascular-flow"}-project-20261022'
 if project=='brava':
  title='BraVa CFD 与微泡项目';entry='computation_examples/brava_flow_roi_18mlmin'
  details=f'''本分支仅保存 BraVa 项目与必要共享依赖。小鼠血管项目见 `sync/mouse-vascular-flow-project-20261022`。两个项目不合并为同一个计算模型；祖先提交仍保留原仓库历史。

当前科学状态以 [{entry}/ACTIVE_DEVELOPMENT_SCOPE.md]({entry}/ACTIVE_DEVELOPMENT_SCOPE.md) 为准：P2/P1 路线在资源门槛处停止；Gmsh V2 场局部精度仍未通过，未运行新正式群体。

- [正式流场身份]({entry}/reports/ACTIVE_FLOW.json)：18 mL/min 历史 balanced_pressure_final；候选 0/15 为同一解的打印姿态。
- [旧正式群体账本]({entry}/population_v1/runs/realization_0000_contactfix/SUMMARY.json)：1565 条接受轨迹，458/527/437 出流，141 观察截止、2 数值失败。原始 states/chunk_*.npz 全部保留。
- [最新四视角]({entry}/population_v1/candidate_0/OPEN_RESULTS.html)与[流场可视化]({entry}/visualization/candidate_0/OPEN_RESULTS.html)。
- [P2 线性代数报告]({entry}/population_v1/analysis/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)。
- [新服务器迁移与独立环境验证](deployment/newserver_20261009/MIGRATION_AND_VALIDATION_ZH.md)；当前入口和激活脚本在 `deployment/current_server/`，环境验证不等于局部物理精度认证。

当前固定灌注协议尚无新正式大批量结果；未来新边界按中心线阻力估算出口压力，历史校准场身份保持不变。服务器当前项目与本地已有对应文件 SHA 全部一致；服务器独有、与当前诊断相关的源码与执行证据从已核验备份补齐。
'''
 else:
  title='小鼠血管流场与配套微泡项目'
  details='''本分支保存小鼠血管项目；BraVa 项目见 `sync/brava-flow-project-20261022`。共同粒子/求解器代码是共享依赖，不代表两个几何或结果相同。祖先提交仍保留原仓库历史。

- [ROI-only 算例](computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/CASE_README.md)：1500 条、dt 0.5 ms，1315 出流（787/61/467）、185 接触支持静止、0 执行失败。
- [best-feasible 算例](computation_examples/mean-2p0-mmps-A-best-feasible-balance-v1/CASE_README.md)：1500 条、dt 0.5 ms，1315 出流（0/133/1182）、185 接触支持静止、0 执行失败。
- [Network-H0 基线](computation_examples/mean-2p0-mmps-A-H0-pressure-v1/CASE_README.md)：500 条、dt 1 ms，430 出流（0/139/291）、70 接触支持静止。
- [原 2 mm/s 基线](computation_examples/mean-2p0-mmps/CASE_README.md)与[早期 FEM Stage-Q 参考](computation_examples/fem_reference_stage_q/CASE_README.md)。
- [WSS 验证报告](computation_examples/_shared/wss_validation_v2/WSS_VALIDATION_REPORT_V2.md)及 `validation_cases/`；细档血管和 O2 压力敏感性继续为 SKIPPED_BY_USER。
- [ROI-only 流场图表](computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/visualization/OPEN_RESULTS.html)和[微泡图表](computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/microbubble/OPEN_RESULTS.html)。

四套冻结场、求解配置与日志、WSS 表面、原始轨迹及验收记录均保留。`geometry_inputs/` 只补入此流场来源使用的源样本、ROI/表面和共享几何代码，不复制整个上游打印模型库。小鼠项目未迁入当前 BraVa 服务器；本次同步依据本地保留的权威结果，没有声称重新核验旧服务器。
'''
 text=f'''# {title}（2026-10-22 同步）

{details}
本次是现有文件同步，不启动或恢复 CFD、轨迹、网格加密或敏感性计算，不修改科学数据。源目录、原有兼容链接和完整本地备份保留。

纳入 {len(inc):,} 个源文件，约 {summary['source_bytes']/1e9:.2f} GB。保留代码、关键日志、冻结网格/流场、原始正式状态、诊断表格和成品图表；排除可安装环境/二进制、缓存、重复中间帧、废弃批次大体积状态、可由已保留输入重建的全域审计数组。排除是发布筛选，没有删除原始数据；普通 Git 不上传单文件 ≥100 MiB 的派生大文件。没有 LFS，也没有将文件切片规避限制。

逐文件来源/大小/SHA256 见 [INCLUDED.json](sync_20261022/INCLUDED.json)；[EXCLUDED.json](sync_20261022/EXCLUDED.json) 记录排除原因、原位置及大文件哈希；[LINKS.json](sync_20261022/LINKS.json) 记录相对兼容链接。少数上游 `.gitattributes` 的 LFS 规则只在发布副本中停用，原文与前后哈希已保留。原生第三方测试库中的既有 LFS 指针仅为上游测试资源引用，不是本项目冻结 CFD 文件。

历史配置和报告中的绝对路径保持原始字节；这不是空白机器一键运行包。HTML/视频可下载本分支后打开；GitHub 页面直接预览部分格式受限。克隆时可加 `--single-branch --depth 1` 仅获取此分支当前快照。

同步分支：`{branch}`。最终提交和远端核验回执在本地 `github_sync_audits/flow_projects_20261022/`。
'''
 (root/'README.md').write_text(text)
 # Source identity and publication safety check. No scientific solver or integration.
 errors=[];findings=[];lfs=[]
 for row in inc:
  p=root/row['path']
  if sha(p)!=row['sha256']:errors.append({'path':row['path'],'error':'hash_mismatch'})
  if p.stat().st_size>=100*1024*1024:errors.append({'path':row['path'],'error':'large_git_blob'})
 for p in root.rglob('*'):
  if not p.is_file() or p.is_symlink() or p.name=='.git':continue
  with p.open('rb') as f:
   head=f.read(8192)
   if head.startswith(b'version https://git-lfs.github.com/spec/v1'):lfs.append(str(p.relative_to(root)))
   if b'\0' in head:continue
   data=head+f.read()
  for rule,pattern in patterns.items():
   if pattern.search(data):findings.append({'path':str(p.relative_to(root)),'rule':rule})
 for row in links:
  p=root/row['path']
  if not p.exists():errors.append({'path':row['path'],'error':'dangling_symlink'})
 forbidden='computation_examples/mean-' if project=='brava' else 'computation_examples/brava'
 if any(x['path'].startswith(forbidden) for x in inc):errors.append({'error':'project_scope_mixing'})
 result={'project':project,'source_hashes_verified':len(inc),'errors':errors,'credential_candidates':findings,'upstream_existing_lfs_pointer_paths':lfs,'status':'PASS' if not errors and not findings else 'REVIEW_REQUIRED'}
 (AUDIT/f'{project.upper()}_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
 (meta/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
 print(project,summary,'validation',result['status'],'errors',len(errors),'secret_candidates',len(findings),'upstream_lfs_refs',len(lfs),flush=True)
