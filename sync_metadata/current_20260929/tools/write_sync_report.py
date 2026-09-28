from pathlib import Path
import json
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260929'
collection=json.loads((M/'collection_summary.json').read_text());server=json.loads((M/'server_summary.json').read_text())
v=json.loads((M/'validation/snapshot_validation.json').read_text());assert v['PASS']
(M/'EXCLUSIONS.md').write_text('''# 排除项与恢复方式

本次保留正式网格、校准/正式CFD的原始结果与日志、失败记录、两姿态10段流场MP4、4K PNG/PDF、微泡样本与所有已生成的基准/未完成轨迹。既有小鼠1500条、RBC及WSS结果不删除。

未纳入的内容：

- Python/Conda环境、Git数据库、安装包、缓存、套接字/PID、可重建共享库/对象文件/求解器二进制。求解器完整修复源码、补丁、构建命令和二进制SHA256保留。
- BraVa `microbubble/data/gpu_mesh_input.npz`，109,638,745字节。它是完整网格/节点场/梯度/壁面数组的派生输入，超过GitHub单文件限制；正式VTU、NPZ场和壁面数据保留，CUDA计算记录也保留。原副本仍在WSL及服务器。
- 以前已经排除的用户取消细档大网格、第三方示例输出和上游未下载的LFS真实载荷。LFS指针文本记录于`server_inventory.json`，不能视为真实测试数据。
- 与本轮流场无关的大型TopBrain原始影像、全数据集、切片安装包、G-code/3MF/ZIP和重复完整盒体；保留所用BG001 SWC、refined ROI、BALANCED候选、四端口芯、姿态清单、配置与开发源码。其它打印候选及打印专用完整输出仍在原`vascular_printing/`工作区。

逐项排除路径、原因、大小和可用SHA256见 `excluded_files.json`、`server_inventory.json`。目录级运行环境不会逐文件扫描私有配置。它们是排除清单，不表示源文件已删除。

仅恢复现有GPU输入文件（不会启动计算）：

```bash
scp -p vast4090:/workspace/brava_flow_roi_18mlmin_20260928/microbubble/data/gpu_mesh_input.npz brava_flow_roi_18mlmin/microbubble/data/
sha256sum brava_flow_roi_18mlmin/microbubble/data/gpu_mesh_input.npz
```

校验值请与清单对应记录核对。不要通过运行整个`campaign.py prepare`恢复缓存，以免新建或覆盖样本契约。当前微泡取消标记应保留，用户明确授权恢复之前不执行计算命令。

同步源/配置/结果保持原始字节；导航文档变换和Git规则归档另有记录。使用者仍需按原路径契约配置科学环境；当前快照不是一键部署包。
''')
(M/'SYNC_REPORT_ZH.md').write_text(f'''# 2026-09-29 当前开发内容增量同步

目标仓库：`OZAKI39/ulm-3d-vascular-model-generation`。按用户指定，更新现有分支 `sync/roi-flow-mb1500-wss-stop-analysis-20260928`，父提交 `ca0ae424d01717f6a231bbea91405cbea8f095a8`。工作副本为 `/home/lzy/projects/temp_storage/github_sync_20260928_latest`；原科学工作区不提交、不改写。

## 新增及更新

1. BraVa candidate 0/15所用四端口芯、真实网格、18 mL/min校准与正式压力出口CFD、全部本轮小算例/失败日志、稳态记录、独立WSS后处理与内部通量诊断。
2. 两姿态各5段流场动画和4K图，实际运行数据及媒体哈希。固定Z轴旋转、配色和标签沿用原风格，显示单位/色标按新数据。
3. 修复GPU主机视图同步的独立求解器源码、补丁、构建记录和CPU/GPU两步对照；950项上游LFS占位单独记录，不冒充真实载荷。
4. 1500个微泡出生样本、CUDA梯度核验、完整单条原/编译核基准及审计、8条未完成生产轨迹、取消记录、未完成缓存优化草稿。**正式完成0条；本批没有微泡动画、没有1500条交付通过标记。**
5. 当前共享几何/ROI源码、配置、指南及所用BraVa源几何；保留原小鼠管网、FEM、1500条轨迹、停止分析、WSS审计和独立RBC内容。

流场与微泡不能混淆：BraVa正式流场SHA256=`{v['flow_sha256']}`；原小鼠1500条仍绑定其自己的流场。当前BraVa内部截面流量存在约百分之几的局部守恒缺陷，独立散度核验也已随附；不因分流接近均匀而隐藏这一限制。

## 来源与同步范围

本地当前来源核对 {collection['source_files']:,} 个常规文件、{collection['logical_bytes']:,} 字节及{collection['symlinks']}项链接记录。逻辑体积包含已在分支中的重复数据，不是本次网络上传量。原网络/小鼠CFD实际开发来源仍采用`temp_storage/github_sync_20260927`，避免用另一处旧Git HEAD覆盖新开发。其余来源及精确映射见`source_files.csv`。

服务器本轮BraVa目录读取{server['inventory_files']}项；文件映射计数为 `{server['actions']}`。修复源码、执行配置和服务器新增报告补入快照，已有相同字节直接映射。未对服务器执行启动/停止求解、环境安装或清理。

同步保留服务器微泡 `STOPPED` 状态；配置 `autostart=false`、`autorestart=false`。`USER_STOP.json`、停止后的主入口保护和未完成记录均保留。代码中存在恢复/渲染/报告函数不表示它们已执行。

## 此次实际核验

| 检查 | 实际结果 | 范围 |
|---|---|---|
| 本地源与同步副本SHA256 | {v['source_files_checked']:,}项；源变更{len(v['changed_sources'])}项 | 文档导航转换单独记录；未改写科学数据 |
| 服务器文件映射SHA256 | {v['server_mappings_checked']:,}项 | 只读远端清单；无新求解 |
| 受保护微泡科学源码 | {v['protected_scientific_sources_verified']}项一致 | 原动力学源文件字节保持不变 |
| 已完成单条基准文件哈希 | {v['pilot_files_hash_verified']}项一致 | 原/编译核轨迹与审计一致；非1500条批次通过 |
| 正式CFD输入、冻结场、WSS及两姿态媒体 | 通过，{v['flow_movies_hash_verified']}段媒体哈希一致 | 原解码记录已保留；本次未再渲染或积分 |
| 微泡取消状态 | `CANCELLED_BY_USER`；0完成、8未完成 | 缓存草稿验证取消，不能写PASS |
| 新任务脚本语法 | {v['syntax_only_files']}份可解析 | 未导入或执行科学流程 |
| 常见凭据模式及单文件体积 | 未发现待处理匹配；无≥100 MiB新文件 | 文件大小和模式检查不等同全面安全认证 |
| 新CFD / 新轨迹积分 | 0 / 0 | 本轮纯同步 |

源路径清单、Git规则归档、导航变换和原文保留于本目录。`snapshot_manifest.csv`与Git索引校验记录由发布前步骤生成；远端分支及提交树在推送后再次核对。旧测试/验收记录保留原日期，本次没有把上轮测试称为新执行。

## 排除与使用

详见[EXCLUSIONS.md](EXCLUSIONS.md)。正式结果及动画保留，主要排除运行环境、二进制、超过单文件限制的派生GPU输入、已取消细档大网格和独立打印/影像大包。排除文件不从WSL或服务器删除。

根[README](../../README.md)与[BraVa状态说明](../../brava_flow_roi_18mlmin/CURRENT_RESULTS_ZH.md)提供入口。此快照仍包含实际WSL/服务器路径契约，不能宣称在空白机器上直接运行；本次没有修改受保护源码来“适配路径”。
''')
print('Sync report and exclusions written')
