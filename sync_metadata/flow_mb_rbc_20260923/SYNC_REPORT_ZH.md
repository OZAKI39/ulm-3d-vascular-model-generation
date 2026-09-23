# 完整血流 / 微泡 / RBC 工作流 GitHub 同步报告

目标仓库：`OZAKI39/ulm-3d-vascular-model-generation`。新分支：`sync/flow-microbubble-rbc-workflow-20260923`。同步日期：2026-09-23。

## 完成范围

本快照整合当前 Particle 工作树、新 2.0 mm/s 流场工作树、原 FEM 工程的旋转可视化与历史证据、早期 FEM、冻结 SonoVue 分布和来源迁移记录。保留所有对应模块源码、永久测试、脚本、配置、必要网格/输入、正式轨迹、最新图像/动画、CSV/JSON/Excel、中文报告和服务器日志/计算来源。

没有更改求解或运动方程，没有重新进行全量 CFD 或批量轨迹积分。原工作目录中的文件未删除；GitHub 主分支不在本次写入范围。

## 分支策略与来源

以已发布 FEM 分支的 `c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2` 为父提交建立独立工作树，再同步当前文件。没有继承本地未发布 Particle 历史中的两个超 100 MiB CSV 对象，因此新分支不会因为这些历史大对象而超出 GitHub 单文件限制。

来源 HEAD、分支和工作区状态见 `source_git_provenance.json`。同步包括尚未提交的相关源文件修改。已有 FEM 冻结依赖闭包完整保留；基线中的数值文件未改写。`source_README.md` 和 `source_workspace_index.md` 保存原入口说明；新入口为仓库根 `README.md`。

## 精简范围

逐文件选择记录共 29,743 条：包含 16,637 条，已由发布基线覆盖 2,355 条，排除 10,751 条。来源和目标的完整 SHA256 比对已通过 16,637 项。单次清点中排除约 6.29 GiB 文件；此数字不包含整体跳过的虚拟环境、构建目录、Git 内部对象和完整原始脑数据集。

排除大型准入缓存、重复轨迹 CSV、逐状态派生升力数组、非必要中间 CFD 状态、多余视频帧、旧可视化和部署压缩包。正式轨迹 NPZ/JSON、最终场与原生检查点、稳态门槛对应的 result_070.vtu、最终 result_071.vtu、最新旋转可视化的标签/验证数据均保留。没有增加 Git LFS 依赖。

详见 `file_inventory.csv`、`snapshot_summary.json` 和 `REBUILD_OMITTED_ZH.md`。历史审计清单可能列出被排除的文件，当前精简快照必须用本次 `workflow_files.sha256` 检查。

## 校验

- 16,637 个所选源文件与复制结果 SHA256 一致；已发布基线没有缺失文件。
- 扫描 12,477 个文本文件，未匹配到私钥、GitHub/AWS/API 令牌或需人工复核的明文凭据赋值。该扫描不是对任意格式敏感信息的完备证明。
- 流场输入、边界条件、实际新场、守恒/稳态证据、WSS、残差语义、RBC/微泡运动：100 项通过。
- 独立升力公式、零值保护、量纲、正式方程未引入升力、近壁模型一致性及图像检查：13 项通过。
- 共 113 项通过。RBC interaction 的已有零时长接触路径产生 1 条 `RuntimeWarning`；未导致断言失败，本次为同步任务，未改写对应模型。
- 3 项历史只读测试使用原机器绝对路径，因此未作为可迁移测试运行；当前快照的源文件/轨迹完整性由本次 SHA256 清单检验。没有声称重跑全部历史回归测试。
- 初次测试集合缺少 `sv_validation` 导入路径，已通过显式设置本仓库 `PYTHONPATH` 修正运行环境；代码未因该问题修改。`run_checks.py` 封装了最终命令。
- 当前最大已暂存文件约 23.16 MiB，所有对象均低于 GitHub 的 100 MiB 单文件上限。新增文件采用 25 MiB 软上限；已有冻结闭包按原样保留。

测试日志为 `pytest_workflow.log`、`pytest_audit.log`；机器结果为 `preflight.json`、`validation_summary.json`。解释器为 Python 3.13.11，复用已有科学环境。

## 模型边界与使用

旧冻结场与 2.0 mm/s 新场分开保留；现有正式微泡轨迹不能解释为新场轨迹。最新正常形态 RBC/MB 同流演示使用理想化 Poiseuille 场；几何形变代理不是完整膜力学模型。剪切升力量级审核没有写回正式运动方程。原阶段的模型范围、有效性限制、时间步/网格收敛缺口和 GPU 等效性结论均未改变。

下载后先运行 `python3 sync_metadata/flow_mb_rbc_20260923/verify_snapshot.py`。安装科学依赖后可加 `--scientific-inputs` 验证冻结 FEM 与完整 SonoVue 契约，再运行 `run_checks.py`。旧报告中的服务器/WSL 绝对路径是历史来源记录，重跑前需按配置迁移。完整目录入口见仓库根 README。

历史文件的全量 `git diff --check` 检出原始 CSV 的 CRLF 与历史格式空白；为保持冻结哈希，未进行格式化。新增入口文档和打包工具另行进行空白检查。
