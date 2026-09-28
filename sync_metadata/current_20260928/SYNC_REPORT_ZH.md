# 最新 ROI 流场、微泡与 WSS 工作流同步报告

日期：2026-09-28。目标：`OZAKI39/ulm-3d-vascular-model-generation`。

新分支：`sync/roi-flow-mb1500-wss-stop-analysis-20260928`；父提交为 `69c6af94eabb87ba6db0823134c4b92b30a069db`。建立独立工作树收集实际文件，不切换或提交原开发工作树，不改变远端既有分支。

## 纳入的内容

1. 当前血管A/ROI/1D0D、ROI-only边界设计，以及新边界下真实3D FEM的源码、网格、配置、日志、原始结果、独立验证与GPU证据。
2. 当前微泡有限尺寸计算模块、1500条原始轨迹及其9000个结果文件、身份清单、0.5ms导出、主动画及3个附加视角。
3. 当前 WSS、原始面片WSS、应变/剪切率、J1分步WSS动画；原WSS审计、第二轮两档敏感性结果及会议疑问位置说明图。
4. 最新中途停止原因分析，包括185个位置的静态复算证据、两张PNG/PDF、中文报告、CSV及全部诊断脚本。
5. 当前保留的H0/best-feasible参考、粒径分布、独立RBC源码/轨迹/展示、HemoCell资料及服务器路径索引。

本地收集 **37,758个常规文件，6,687,884,208字节**；另保留90个相对链接及23个原Git规则文件。这是逻辑文件总量，包含重复输入/参考，Git按内容去重后上传量不同。最终Git有效文件及体积以 `validation/snapshot_validation.json` 和 `snapshot_manifest.csv` 为准。

## 来源优先级

- 最新管网及ROI CFD实际开发在 `/home/lzy/projects/temp_storage/github_sync_20260927/` 下；这一目录已经是活跃工作树，不能只从 `/home/lzy/projects/ulm_3D_vascular/` 取旧版本。此次保留其全部新增代码、报告和算例。
- 最新微泡来自 `/home/lzy/projects/ulm_particle_formal_p9a5/`。
- 网格工具、WSS审计和最新动画来自 `/home/lzy/projects/formal_3D_flow_solver/`。
- 几何共享代码按前次整理范围从当前 `vascular_printing/` 同名文件刷新；源目录/仓库目录逐文件映射见 `local_file_inventory.csv`。

发现一处同路径不同版本：`ulm_3D_vascular/utils/cfd_preprocess/one_d_flow.py`。采用当前管网开发工作树的实现，另一处旧副本及两者哈希记录于 `source_precedence.json`。未拼接或改写科学实现。

唯一对收集文件进行的内容转换是根 `ACTIVE_VASCULAR_WORKFLOW.md` 中指向旧工作树的目录导航链接，改为当前仓库内路径；原文本另存，变换与双哈希见 `documentation_transforms.json`。符号链接转换为快照内相对路径；科学源码/配置/结果保持原字节。原Git规则归档，根 `.gitattributes` 关闭LFS/EOL变换。

## 服务器只读核对

通过已有 `vast4090` 别名核对当前ROI流场/微泡、上一批流场/微泡、RBC源码与结果、求解器源码及构建记录。**27,005个文件，2,911,634,257字节**纳入远端哈希清单；未缺失任何指定根目录。

按内容去重后只补入 **6个文件，656,478字节**，包括上一批服务日志和预览媒体；最新ROI生产文件均已在本地快照中找到相同字节。950个上游LFS指针映射到明确标注的指针元数据，其余映射均有实际文件。完整映射见 `server_file_map.json`，服务器清单和求解器二进制身份见 `server_inventory.json`。

未发送启动/停止求解指令、未更改服务器配置或数据。

## 本次实际运行的检查

| 检查 | 结果/证据 | 边界 |
|---|---|---|
| ROI-only与参数化管网/边界交接测试 | 106通过，0失败；`validation/network_tests.log` | 720条VTK/NumPy弃用警告，未导致失败；没有启动CFD |
| 原文件与同步副本逐字节核对 | `validation/snapshot_validation.json` | 包括文档变换、相对链接与源文件未变检查 |
| 最新微泡交付/逐轨迹哈希 | 同上，1500条、9000个结果文件及114个受保护源码 | 验证快照内文件；不是重新积分 |
| 最新流场与微泡身份一致 | VTU SHA256=`fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12` | 不代表再做一次CFD收敛实验 |
| 服务器文件映射 | `server_summary.json`及快照验证 | 仅哈希、只读收集 |
| 发布内容检查 | 文件大小、相对链接、常见密钥模式、Git索引字节核对 | 索引检查另存 `validation/staged_validation.json` |

现有数值验收、WSS数值检查、动画检查及停止分析均随结果原样保留；本轮没有把历史验收写成新运行。细档真实血管CFD与出口压力敏感性保持用户取消状态，不追加求解。

## 复核及复现

`snapshot_manifest.csv` 记录文件/链接的大小、SHA256和Git blob ID；自身及动态检查输出不自引用。`tools/validate_snapshot.py` 需要原本地目录以检查源文件未变。纯下载副本可用清单逐文件计算SHA256；链接按链接文本计算哈希。

本次测试命令（从仓库的 `ulm_3D_vascular/` 运行）：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -B \
  -m pytest -q -p no:cacheprovider tests/network_roi_only tests/network_parameterized
```

收集脚本位于 `tools/package_local.py`、`tools/collect_server.py`。它们只应在独立同步副本内执行，依赖原WSL路径和SSH访问；不会启动科学计算。原始微泡和各可视化目录中的复现说明保持其已有路径契约，尚不能宣称在任意空白机器上直接运行。运行环境、第三方LFS真实载荷及排除的大文件不在本分支中，见 [EXCLUSIONS.md](EXCLUSIONS.md)。

根 [README](../../README.md) 提供最新入口；上一次同步报告保留在 `sync_metadata/current_20260927/`，其数值和状态仅描述当时版本。
