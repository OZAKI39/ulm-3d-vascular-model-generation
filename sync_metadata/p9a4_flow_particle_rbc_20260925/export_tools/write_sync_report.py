from pathlib import Path
import json,shutil
A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination']);M=D/'sync_metadata'/C['snapshot']
text='''# P9-A.4 / 3D FEM / Microbubble / RBC 同步记录

本分支是当前本地与服务器开发资料的发布快照。用户已明确授权同步到GitHub新分支；这里的发布发生在P9-A.4开发验收之后，因此原阶段报告中的“没有push”保留为当时的历史记录。

- 仓库：`OZAKI39/ulm-3d-vascular-model-generation`
- 新分支：`sync/p9a4-flow-particle-rbc-20260925`
- 基线：`sync/network-h0-particle-cleanup-20260925`，`abb3ab05fb8dcddcf120765702b682f23a63a71d`
- 采集开始：2026-09-25 00:17:43 UTC；这是一段时间内只读采集形成的快照，不宣称对所有正在开发的机器执行了原子冻结。
- 发布工作树：`/home/lzy/projects/github_sync/p9a4_flow_particle_rbc_20260925`
- 本机同步审计/传输包：`/home/lzy/archives/github_sync_p9a4_20260925T001743Z`
- 服务器采集归档：`vast4090:/workspace/archives/github_sync_p9a4_20260925T001743Z`

## 保存的开发内容

| 类别 | 内容与入口 |
| --- | --- |
| 血管网络与ROI | [vascular_network](../../vascular_network/)，包括A-Network H0、ROI边界、网络压力/流量分配及当前第一方源码 |
| 3D FEM | [formal_3D_flow_solver/FEM_SimVascular](../../formal_3D_flow_solver/FEM_SimVascular/)，包含算例、网格、入口/出口配置、固定流场、求解脚本、报告与既有日志 |
| 微泡与RBC | [particle_3d](../../particle_3d/)，保留微泡动力学/轨迹代码、RBC运动/变形代理、测试、模型合同及历史验收资料 |
| 最新P9-A.4 | [particle9a4_population_inlet](../../particle_3d/reports/particle9a4_population_inlet/)，321个新增代码/合同/测试/报告/图件/日志/数据文件逐字节保存；3份大文本以gzip表示 |
| 输入分布与可视化 | SonoVue固定分布、旋转流线/速度/压力/WSS可视化、30泡输出及8组图件 |
| 服务器证据 | [server_evidence](../../server_evidence/)，新增8,350份文件；另29,772条来源记录与已有文件内容完全一致，通过SHA映射定位 |
| 路径索引 | [WORKFLOW_PATHS_ZH.md](../../WORKFLOW_PATHS_ZH.md)及[path_inventory](path_inventory/)，保存本机/服务器原始路径和目录扫描 |

既有快照中46,671个受版本控制的文件为基础。本次重新核对29,159份原始已选文件：29,149份相同、9份当前TopBrain/BRAVA源码已变化、1份旧模块在原工作区已被移除。随后单独扫描479份当前vascular第一方源码/配置/测试/文档，又采集14份新增/更新文件。已从当前导出目录移除原工作区不再存在的`partial_fgw_transfer.py`；其历史版本仍在父提交中。详见[inherited_source_recheck.json](inherited_source_recheck.json)、[vascular_source_refresh.json](vascular_source_refresh.json)。

TopBrain/BRAVA是独立开发线，只保存当前源码快照；本次165项FEM/Particle验证不代表该独立开发线或全部历史HemoCell/LAMMPS代码已经重新验收。其完整上游数据集未重复上传。

## 当前科学状态与输入角色

NEW H0流场SHA256：`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。本次没有重算CFD或生成新轨迹；同步已有经过审核的结果。Particle默认frozen_reference仍属于历史OLD equal-pressure流场，NEW必须依照P9-A.4/Network-H0合同显式选择。

P9-A.4保留100,000个source proposals、19,221个accepted births和80,779个拒绝事件，拒绝记录未删去。固定浓度平台是明确的建模假设，并非实测绝对浓度。NEW H0首30泡状态为28 completed、2 supported stationary、0 solver fail。当前阶段为`P9A4_POPULATION_INLET_READY_FOR_LARGE_SAMPLE_REVIEW`，没有运行正式500。原P9-A.1科学接口及旧入口合同保持原始字节。

历史RBC Poiseuille展示、HemoCell与LAMMPS资料保留各自来源，它们不是NEW H0流场上的RBC生产轨迹。Taylor-Hood验证保留为已经停止的历史阶段。

## 大文件处理与服务器去重

重要科学大文本没有舍弃。以下3份原始文件用确定性gzip无损保存，记录原始/压缩SHA与大小于[packed_artifacts.json](packed_artifacts.json)：

| 原始文件 | 原始字节 | 压缩字节 |
| --- | ---: | ---: |
| `data/inlet100k/accepted_births.json` | 19,524,555 | 3,045,855 |
| `data/inlet100k/proposal_ledger.jsonl` | 92,108,810 | 18,192,341 |
| `data/legacy_b_inlet_500.json` | 61,236,773 | 4,119,054 |

以上相对路径位于`particle_3d/reports/particle9a4_population_inlet/`。运行下述命令可恢复旧报告和checkpoint读取时使用的原路径，不需修改科学代码：

```bash
python3 sync_metadata/p9a4_flow_particle_rbc_20260925/restore_large_artifacts.py
```

还原工具同时检查gzip和原始内容的SHA；遇到不同内容的既有文件拒绝覆盖。同SHA文件允许重复执行。还原出的3份原始文件被精确列入.gitignore，避免重复提交。已在临时目录实际还原、用原始`PopulationLedger.restore`恢复checkpoint并逐项核对全部proposals/births；见[packed_restore_validation.json](packed_restore_validation.json)。

服务器新增实际文件约920.81 MiB。29,772条相同内容来源通过`representation.path`引用已有文件；引用若是gzip，则明确标记`GZIP_OF_IDENTICAL_BYTES`。该映射经过解码后的SHA核对。5个历史缺失LFS对象的指针只保存为JSON证据，不伪装成已拿到实体数据。继承快照中的历史LFS说明保持原状；本次仓库文件使用普通Git对象，不新增LFS依赖。

[server_delta_inventory.json.gz](server_delta_inventory.json.gz)列出每份已扫描服务器文件的原路径、状态、大小、SHA/排除原因；[server_delta_summary.json](server_delta_summary.json)给出统计。[server_reference_validation.json](server_reference_validation.json)记录对实际发布目标的最终检查。`EXCLUDED_FROM_DELTA`表示不从此次服务器增量再次纳入，不表示删除继承快照中已经保存的资料。

排除环境、安装包、编译结果、缓存、完整PETSc/Palabos/TBB等未修改第三方安装/源码树、完整上游脑数据、大型渲染帧及非必要中间/重复状态。服务器传输包检查后剔除了15,335份这类第三方文件；原服务器未删除任何文件。必要第一方代码、已有固定第三方版本/补丁、环境版本记录及运行证据保留。

## 验证与复现

在同步工作树内执行portable Particle/Network/Flow 115项、P9-A.4 35项及legacy入口/open-cap 15项，共**165项通过，0失败、0错误、0跳过**。最后一次运行64.37秒。完整日志、JUnit和摘要见[validation](validation/)；依赖版本沿用[requirements-verified.txt](../network_h0_particle_20260925/requirements-verified.txt)。

```bash
python3 sync_metadata/p9a4_flow_particle_rbc_20260925/verify_snapshot.py --scientific-inputs
python sync_metadata/p9a4_flow_particle_rbc_20260925/run_checks.py --output /tmp/ulm-p9a4-checks-new
```

测试输出目录必须是全新目录。校验器检查本分支所有发布文件的SHA、3份压缩原始文件，以及沿用的40份网络输入、114份受保护科学源码、114份server bundle源码和OLD/NEW冻结流场合同。`SNAPSHOT_SHA256.txt`不包含自身。旧快照的整仓manifest仍是历史记录，新的全仓内容请使用本目录的校验器。

321份P9-A.4文件和源工作区进行字节级比对的结果见[p9a4_byte_validation.json](p9a4_byte_validation.json)。发布前扫描结果见[preflight_audit.json](preflight_audit.json)；其中安全扫描针对新增/修改内容中常见凭据格式，不能解释为对所有可能凭据格式的绝对证明。

## 工作区与文档处理

发布只在独立worktree上创建提交。原开发工作区的HEAD、branch、index及status前后记录见[source_git_before.json](source_git_before.json)、[source_git_after.json](source_git_after.json)；独立开发文件若在采集后继续变化，记录于[concurrent_source_changes.json](concurrent_source_changes.json)，不回滚、不覆盖。服务器只读采集，另写入新归档目录。

本次更改README、导航、忽略规则与同步工具，使新入口及gzip还原方式可见；原文副本在[original_docs](original_docs/)，修改列表在[packaging_edits.json](packaging_edits.json)。原科学报告、审计合同和运行日志保持原文，包括其历史机器绝对路径、旧状态和当时的“未push”说明。export_tools保存本次打包过程的来源脚本，依赖采集归档中的context，不作为独立科学复现入口。

新分支独立发布，不合并main，不提交原开发工作区的暂存内容。提交SHA及远端验证收据保存在本机归档，避免在manifest内自引用当前提交SHA。
'''
(M/'SYNC_REPORT_ZH.md').write_text(text)
shutil.copy2(A/'context.json',M/'snapshot_context.json')
shutil.copy2(A/'audit/local_summary.json',M/'local_collection_summary.json')
for p in (A/'scripts').glob('*.py'):
 if p.name!='server_delta_payload.py':shutil.copy2(p,M/'export_tools'/p.name)
print('SYNC_REPORT_READY')
