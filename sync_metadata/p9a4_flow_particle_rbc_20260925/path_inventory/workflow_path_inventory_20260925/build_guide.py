from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import json,csv,re,subprocess
W=Path('/home/lzy/projects');OUT=W/'workflow_path_inventory_20260925'
keys=json.loads((OUT/'key_paths.json').read_text());remote=json.loads((OUT/'remote_paths.json').read_text());local=json.loads((OUT/'local_paths.json').read_text())
def add(loc,cat,why,path):
 p=str(path)
 if loc=='WSL':exists=Path(p).exists();kind='directory' if Path(p).is_dir() else 'file' if Path(p).is_file() else 'missing'
 else:
  exists=p in rdirs or p in remote['rbc_related_paths'];kind='directory' if p in rdirs else 'file'
 keys.append(dict(location=loc,category=cat,purpose=why,path=p,exists=exists,kind=kind))
rdirs={x['path'] for x in remote['directories']};srv='root@50.115.148.16:4159'
A=W/'ulm_particle_population_inlet_p9a4';P=W/'ulm_particle_3d_particle0/particle_3d'
add('WSL','历史粒子使用的冻结 FEM','当前默认frozen_reference内的历史equal-pressure场（OLD SHA）',P.parent/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference/flow/steady_flow_mean_2p0_mmps.vtu')
add('WSL','Git 同步快照','共享Git对象库（各worktree依赖）',W/'ulm-3d-vascular-model-generation/.git')
for name in ['vascular_network','particle_3d','formal_3D_flow_solver','sonovue_size_distribution_v0','server_evidence','source_dependencies','THESIS_SUMMARY']:
 add('WSL','当前 P9-A.4 集成工作树内容',name,A/name)
H='/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1'
for name in ['SV_MESH','reports','run','run/1-procs']:add(srv,'服务器当前主线','H0算例 '+name,H+'/'+name)
N='/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z'
for name in ['particle_3d/src','inputs','scripts','outputs/NEW','outputs/NEW_POINT','logs','final_review']:add(srv,'服务器当前主线','Network验证 '+name,N+'/'+name)
R='/workspace/particle9a4_population_inlet_20260924T233708Z'
add(srv,'服务器RBC代码副本','当前Particle部署携带的RBC模块；不代表本次执行过RBC轨迹',R+'/particle_3d/src/particle_3d/rbc.py')
for stem in ['hemocell_restore','microbubble_lammps','lammps_active']:
 cat='服务器历史 HemoCell / LAMMPS 独立工程'
 for row in remote['directories']:
  if row['root']=='/workspace/'+stem and row['depth']==1:add(srv,cat,stem+' '+Path(row['path']).name,row['path'])
for name in ['work/rbc_stage1_20260915_165632','results/rbc_stage1_20260915_165632','logs/minimal_migration_20260915_143631']:
 add(srv,'服务器历史 HemoCell / LAMMPS 独立工程','HemoCell RBC '+name,'/workspace/hemocell_restore/'+name)
for name in ['github_sync_network_h0_20260924T221559Z','vascular_workflow_cleanup_20260924T214356Z']:
 add(srv,'服务器归档','同步/退役代码归档','/workspace/archives/'+name)
# Add every report/output stage without trying to infer physics solely from its name.
for group in ['reports','outputs']:
 for q in sorted((P/group).iterdir()):
  if q.is_dir() and not any(k['path']==str(q) for k in keys):add('WSL','Particle 完整阶段目录',group+'/'+q.name,q)
with (OUT/'key_paths.json').open('w') as f:json.dump(keys,f,indent=2,ensure_ascii=False)
with (OUT/'key_paths.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['location','category','purpose','path','exists','kind'],extrasaction='ignore');w.writeheader();w.writerows(keys)
with (OUT/'remote_rbc_paths.txt').open('x') as f:f.write('\n'.join(remote['rbc_related_paths'])+'\n')
order=['当前 P9-A.4','当前 P9-A.4 集成工作树内容','当前 vascular A → Network H0 → FEM','Network-H0 配对验证与冻结动力学','微泡和 RBC 共用计算工程','微泡计算与轨迹关键位置','RBC 源码与验证入口','RBC 运动、形状与展示','服务器当前主线','服务器RBC代码副本','服务器历史运行/独立模拟工程','服务器历史 HemoCell / LAMMPS 独立工程','服务器归档']
groups=defaultdict(list)
for k in keys:groups[k['category']].append(k)
for cat in groups:
 if cat not in order:order.append(cat)
header=f'''# 血管流场、微泡与 RBC 全流程路径索引（2026-09-25）

核查时间 UTC：本地 {local['checked_utc']}；服务器 {remote['checked_utc']}。本地路径逐项存在性检查，服务器通过 SSH 只读扫描，主机 `{remote['hostname']}`，连接 `ssh -p 4159 root@50.115.148.16`（本地别名 `vast4090`）。本次仅新增路径索引和扫描记录，没有运行计算、修改科学代码或移动结果。

**先区分版本**：当前正式场是 NEW Network-H0 outlet-pressure 场，SHA `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。最新入口开发在独立 P9-A.4 worktree，100k source audit与30泡smoke已完成；本轮未推送的新增P9-A.4代码/结果只在新worktree及对应服务器目录，不能假定已包含于GitHub同步分支。P9-A.4状态为READY_FOR_LARGE_SAMPLE_REVIEW，正式500未运行。

原Particle工作树的默认frozen_reference仍是OLD equal-pressure 2 mm/s场（SHA `129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d`）；Network配对runner及新P9-A.4入口显式选择NEW。更早P8.2A/PPT500使用更早Stage Q冻结场，不能与新H0结果混用。当前rotate_visualization的输入manifest已经绑定NEW H0，尽管其中兼容文件名仍写mean_2p0_mmps。

RBC与微泡共用particle_3d，不存在单独的当前RBC顶层仓库。已有RBC–MB flow_rotation演示采用理想化Poiseuille场；不是NEW Network-H0 RBC生产轨迹。rbc_mb_hydrodynamic有用户停止标记。服务器另有历史HemoCell RBC和LAMMPS/Palabos工作目录，下文单列；本次只核查其路径，不替这些历史模型作新的有效性验收。

**清单范围**：{len(keys)}条关键路径；本地目录扫描{len(local['directories'])}条、服务器目录扫描{len(remote['directories'])}条。目录扫描最多向各工程根下展开5层，环境、Git对象、第三方源码/构建和大原始数据内部不作递归展开；不是全磁盘逐文件清单。RBC相关命名额外检查至8层，共{len(remote['rbc_related_paths'])}条，包含多个部署副本。

完整表格：[关键路径CSV](workflow_path_inventory_20260925/key_paths.csv)、[本地目录CSV](workflow_path_inventory_20260925/local_directories.csv)、[服务器目录CSV](workflow_path_inventory_20260925/remote_directories.csv)、[远端RBC相关路径](workflow_path_inventory_20260925/remote_rbc_paths.txt)。需要逐个文件时可查最新同步快照的local_file_inventory.csv / server_file_inventory.json，再加P9-A.4的delivery_manifest.json；这些是各自快照时间的清单，不冒充当前全盘实时清单。

目录中的 `src/`、`scripts/` 是源码/入口；`contracts/` 是模型与输入合同；`data/`、`outputs/`、`run/` 是数据/轨迹/原生结果；`reports/`、`logs/` 保存报告及执行记录；`figures/`、`animations/` 保存展示。具体以各阶段目录为准。
'''
lines=[header]
for cat in order:
 lines+=['',f'## {cat}','','| 用途 | 已核对路径 | 状态 |','| --- | --- | --- |']
 for k in groups[cat]:
  p=k['path'];state='存在' if k['exists'] else '**旧索引路径已失效**'
  link=f'[{p}](<{p}>)' if k['location']=='WSL' and k['exists'] else f'`{p}`'
  lines.append(f'| {k["purpose"]} | {link} | {state} |')
lines+=['','## 已知旧索引差异','',
 '旧索引中的 `/home/lzy/projects/ulm_particle_3d_particle0/formal_3D_flow_solver/FEM_SimVascular/frozen_reference/flow/steady_flow_stage_sv1_3q.vtu` 当前不存在；该默认冻结目录的现存文件为 `steady_flow_mean_2p0_mmps.vtu`，对应OLD equal-pressure场。原Stage Q文件副本与历史结果应按各自manifest查找，不能通过改文件名当作NEW场。',
 '旧WORKFLOW_INDEX_ZH.md / WORKFLOW_PATHS_ZH.md保留原文，其中“最新”只代表其写作日期。旧ACTIVE_VASCULAR_WORKFLOW.md尚未涵盖刚完成的P9-A.4。本次新索引优先说明版本映射，未覆盖旧索引或受保护报告。',
 '退役代码原路径与归档位置见 `/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/audit/local_manifest.json` 和 `server_manifest.json`；代码已退役不表示其历史报告/轨迹已删除。',
 '', '## Git工作树记录','', '```text',(OUT/'git_worktrees.txt').read_text().strip(),'```',
 '', '同一Git对象库被多个worktree引用；这些路径是不同版本/开发现场，不是可以直接互换的重复目录。']
target=W/'WORKFLOW_PATHS_20260925_ZH.md'
with target.open('x') as f:f.write('\n'.join(lines)+'\n')
# Verify all local links and the advertised source/current result paths.
for link in re.findall(r'\]\(<?([^)>]+)>?\)',target.read_text()):
 if '://' not in link:assert (target.parent/link).exists(),link
summary=dict(guide=str(target),key_paths=len(keys),existing_key_paths=sum(x['exists'] for x in keys),missing_paths=[x['path'] for x in keys if not x['exists']],local_directories=len(local['directories']),remote_directories=len(remote['directories']),server_hostname=remote['hostname'],checked_utc=datetime.now(timezone.utc).isoformat(),scope='Read-only filesystem path audit; bounded directory depth; no solver or trajectory execution')
(OUT/'summary.json').open('x').write(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(summary,indent=2,ensure_ascii=False))
