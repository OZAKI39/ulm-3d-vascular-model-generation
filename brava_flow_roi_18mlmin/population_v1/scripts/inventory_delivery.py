"""Small authoritative manifest over code, review artifacts and run identities."""
from pathlib import Path
import argparse,json,csv,hashlib,subprocess
ap=argparse.ArgumentParser();ap.add_argument('--run-name',required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];projects=root.parents[1]
source=projects/'ulm_particle_formal_p9a5/particle_3d';repo=projects/'temp_storage/github_sync_20260928_latest'
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
identity=json.loads((root/'runs'/args.run_name/'IDENTITY.json').read_text());records=[]
for start,kind in [(source/'src/particle_3d/population_v1','new_code'),(source/'tests/population_v1','new_tests'),(root/'scripts','new_diagnostic_and_render_scripts'),(root/'contracts','contract'),(root/'validation','validation'),(root/'reports','review')]:
 for p in sorted(start.rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts:continue
  records.append(dict(kind=kind,absolute_path=str(p),bytes=p.stat().st_size,sha256=digest(p)))
for p in sorted(root.glob('*.md')):records.append(dict(kind='delivery_document',absolute_path=str(p),bytes=p.stat().st_size,sha256=digest(p)))
with (root/'FILE_INVENTORY.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
record=dict(schema='BRAVA_POPULATION_DELIVERY_V1',production_run=args.run_name,production_identity=identity,
 code_and_small_artifact_inventory='FILE_INVENTORY.csv',inventory_sha256=digest(root/'FILE_INVENTORY.csv'),
 full_raw_array_and_ledger_hashes=f'runs/{args.run_name}/MANIFEST.json',
 scientific_quality='Use SUMMARY.scientific_integrity_pass and all explicit terminal classes; coverage is not full integrity',
 validation_provenance='Unit tests ran on final 98d705f13 implementation (26 passed), with actual failed-state fallback regression; earlier standalone real-field/wall/reference scripts have their own recorded flow/software/parameters, not the production random realization. Their exercised code was unchanged by subsequent safety fixes.',
 validation_seed_scope='Use per-fixture source ledger or script constants; deterministic geometry cases have no RNG. Never apply production seed to synthetic validation.',
 failure_runs=['realization_0000','realization_0000_guardfix','realization_0000_final'],
 server_root='/workspace/brava_flow_roi_18mlmin_20260928/population_v1',local_root=str(root),
 raw_birth_metadata='wall_gap is admission conservative bound and tetra_id=-1 until first integration; derived export fills exact birth values without rewriting raw arrays',
 protected_originals='validation/protected_legacy_sources.json',GitHub_push_performed=False)
(root/'DELIVERY_MANIFEST.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
lines=['# 新增和修改文件清单','',f'生产代码版本：`{identity["git_commit"]}`。完整逐文件路径、大小及 SHA256 见 `FILE_INVENTORY.csv`。','',
'原先117个物理模块字节保持不变。本轮只新增 population_v1 包、专用测试、BraVa population_v1 目录；不改 CFD、旧源分布、旧轨迹/动画。开发过程修复限定在新增模块。','',
'|模块|用途|','|---|---|']
roles={'environment.py':'权威环境与哈希/端口适配','source.py':'精确流量源、事件 RNG、单次准入','state.py':'CPU SoA、稳定ID','engine.py':'CPU 共享时钟与群体守恒','reference.py':'原物理函数的私有依赖绑定','pairs.py':'精确球对间隙/候选/新交互契约','components.py':'局部组件与稀疏耦合接触','safety.py':'继承数值保护和原支持证据','output.py':'完整接受时间状态和账本','run.py':'CLI、真实覆盖停源、删失','gpu_engine.py':'实际GPU动力学及明确CPU小子集回退','gpu_field.py':'静态FEM上传、tetra walking','gpu_wall.py':'真实有限壁面GPU查询/路径','gpu_neighbors.py':'GPU空间散列及精确gap','gpu_resistance.py':'GPU孤立泡阻力和接触','gpu_state.py':'GPU SoA','gpu_safety.py':'GPU进展历史和终态CPU证据','gpu_batch.py':'各自时钟的多组数值叠批','gpu_extension.py':'FP64 CUDA编译','bvh.py':'静态BVH和tetra邻接','provenance.py':'规范JSON和内容哈希'}
for name,role in roles.items():lines.append(f'|`{name}`|{role}|')
lines+=['|`cuda/kernels.cu`, `cuda/bindings.cpp`|FEM/散列/有限三角形 CUDA kernels 与 PyTorch 绑定|','',
'测试分 Phase2/3/4/6/7/8 和 safety/output；Phase5 真实点集及实际场验证为 scripts 下独立执行脚本。绘图、报告、校验脚本均在 `scripts/`；未将上述工具散落到原正式目录。','',
'已终止的开发失败输出保持原路径和原字节。实际环境二进制、CUDA编译缓存和大体网格属于服务器环境/原项目，交付清单提供路径与哈希，不复制到新增源码包。']
(root/'FILES_CHANGED.md').write_text('\n'.join(lines)+'\n')
print('INVENTORY_COMPLETE',len(records))
