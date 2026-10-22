# 求解软件安装与旧资料清理报告

日期：2026-10-22。范围：本地独立小鼠血管工程。

## 本次做成了什么

**网格和流场求解软件已独立安装，能够实际完成小规模计算。** 安装根目录统一为 `/home/lzy/projects/new_vessel_flow_solvers`。工程内 `configs/environment.json` 已指向这套安装。此前缺的是这组原生软件，而不是 Python 分析程序；现在这项安装工作已补齐。

**你指定的五类旧资料已从新工程删除，约 2.11 GB。** 当前只提供 ROI-only 正式算例入口，旧方案比较图和链接也已更新。

## 怎样确认软件能工作

| 检查 | 实际结果 |
|---|---|
| 生成一个小网格 | SimVascular 原生 TetGen 成功生成 454 点、2053 个四面体 |
| 小管道流场 | 原 svMultiPhysics 有效源码重新编译；CPU/GPU 分别完成 55 步、5040 单元的小管道计算 |
| 新旧计算对照 | 新速度与旧参考相对差异约 3.2×10⁻¹⁶，压力约 4.4×10⁻¹⁴；均远小于原定 10⁻⁷ 验收限 |
| CPU/GPU 对照 | 速度相对差异 2.73×10⁻¹⁶，压力 1.90×10⁻¹⁶；GPU 运算记录已核验 |
| 小表面处理 | VMTK 圆管法向与封口成功，无开放或非流形边 |
| 小管腔重建 | Ultraliser 非对称分叉小模型成功，最终 3049 点、6094 面，无开放边 |
| 当前方法的小测试 | 23 项通过，0 失败、0 跳过；是下列已有测试集的子集，不重复累计 |
| 工程测试集 | 删除用户指定的三项旧测试后重新运行：114 项通过，0 失败、0 跳过 |
| 最新正式结果 | 13,955 个保留文件的 SHA 全部一致；保存流场、1500 条微泡记录、壁面量及当前源码身份核对通过 |

没有重算完整小鼠血管或大批微泡；符合本次“小规模测试即可”的验收要求。上述结果验证安装和已测功能，不给完整血管的局部物理精度增加新结论。

新软件仍使用原 CFD 方程、原有效补丁和求解选项。系统 CUDA 为 12.6，原服务器记录为 13.2，因此新编译文件拥有新身份；没有冒充原文件。旧服务器求解器/PETSc 身份与两份原环境记录另行保留。运行入口会校验安装清单和文件 SHA，不能任意换程序仍沿用原通过记录。

小测试 CPU 用时 17.29 秒、GPU 61.75 秒。这种极小算例上 GPU 启动和数据搬运开销明显，不代表大算例的速度排序。

## 具体删除与保留

删除的新工程目录：

- `results/best_feasible`：另一套压力方案。
- `results/network_h0`：早期 H0 正式结果。
- `results/mean_2mmps`：原始 2 mm/s 独立算例结果。
- `results/stage_q`：早期 Stage-Q 结果。
- `archive`：精选历史说明和旧执行脚本。

五组数据共 17,465 个文件、2,109,803,234 字节；另移除四个退休算例的入口配置，并更新比较图。删除前逐文件核对 SHA，记录保存在[删除执行凭据](audit/authorized_deletion_receipt.json)和[逐文件清单](audit/authorized_deletion_manifest.json)。

最新 ROI-only 流场、网格、1500 条轨迹、原图表、计算日志、原始血管、必要共享输入和现方法验证依据均保留。原轨迹仍为 O1/O2/O3 出流 787/61/467，185 条接触支持静止、0 执行失败；静止不能记为出流。正式流场 SHA 仍为 `fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12`。

删除前保护了四份很小的必要记录：H0 solver.xml、500 条旧队列的测试输入、两份原求解软件环境记录。它们放到当前测试输入/安装依据目录；保留理由和 SHA 见[保护记录](audit/preserved_before_approved_deletion.json)。名称含 H0/v1/v2/2mmps 的共享输入仍有被当前方法读取的文件，不能误删。

当前数据清单是 `configs/active_data_manifest.json`。原迁移清单、旧 delivery 报告和人工审查文件表保留为历史审计；不会把它们改成“当时只迁入最新结果”。首轮完整性检查发现四份已删除算例说明仍在当前清单中，已更正并重新完整核验，记录见 `audit/retired_case_manifest_correction.json`。

旧工作区、共享文件夹、兼容链接及 BraVa 项目未删除；本次没有访问服务器或推送 GitHub。

## 尚存限制

旧 Ultraliser 对两个特殊极小输入有问题：单段输入在统计函数越界崩溃；完全对称微小分叉在网格后处理时崩溃。后续调试已确认直接原因：相邻面片查找失败后仍读取空结果；异常相邻关系最早在哪一步产生尚未定位，不能简单归因于对称性。证据位于安装根目录 `tests/ultraliser_y_diagnosis/`。普通非对称分叉小测试使用同一程序和原参数成功。失败日志保留，不能把普通小样例通过当作这些旧问题已修复。正式几何未修改。

用户后续明确要求删除的三项旧测试已从当前测试集移除，对应缺失样本入口一并清理。同文件中其他仍可运行的测试保留。删除不记作通过，旧运行日志仍是历史记录；当前测试状态见 `validation/TEST_STATUS.json`。几何工具问题的通俗说明见[当前状态](CURRENT_STATUS_PLAIN_ZH.md)。

## 从哪里查看和使用

- [当前小测试报告](../validation/installed_current_small_check/report.json)
- [当前完整测试日志](audit/after_legacy_test_removal.log)
- [最终数据完整性报告](../validation/installed_stack_integrity_final/validation.json)
- [最新结果图表和原动画入口](../results/current_review/index.html)
- [当前软件配置](../configs/environment.json)
- 安装目录的 `README.md`、`manifests/installation.json`、`manifests/build.json` 和 `tests/`：软件版本、真实运行日志、重建说明。

日常小检查（输出目录每次用新名字）：

```bash
cd /home/lzy/projects/new_mouse_vessel_flow
.venv/bin/python scripts/check_current_workflow.py --output validation/my_new_small_check
```

该命令会执行当前已有小测试、读取保存结果、核验安装身份与原生小测试凭据，不自动重跑 CFD。以后明确启动新 CFD 时，`compute-cfd` 的 `--build-json` 使用 `/home/lzy/projects/new_vessel_flow_solvers/manifests/build.json`；准备和启动命令见工程 README。原生编译脚本与环境均在指定安装根目录。
