# BraVa 当前开发快照同步 — 2026-10-06

目标仓库：`OZAKI39/ulm-3d-vascular-model-generation`。新分支：`sync/brava-p2-algebra-20261006`。

本轮只同步已存在的代码与证据，没有运行CFD、微泡轨迹或新群体。来源为本地BraVa独立仓库 `a4b238b`；从上次已发布BraVa分支的 `d9ece125ef92ac87b87626fdf2a9e474fdb071f2` 建立新分支，将当前工作树文件按清单同步。独立算例Git历史及同步checkout中未发布的旧提交不直接并入；既有本地分支、远端分支和main不改写。

[最新开发报告](../../brava_flow_roi_18mlmin/population_v1/analysis/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md) · [三张300dpi核查图](../../brava_flow_roi_18mlmin/population_v1/analysis/p2_linear_algebra_rerun_20261006/OPEN_RESULTS.html) · [逐文件清单](INCLUDED_FILES.json) · [发布检查](SYNC_VALIDATION.json)

## 包含范围

- 当前BraVa已跟踪源码、adapter、脚本、配置、报告、机器结果、图表、必要NPZ和VTK数据。包括原速度事件定位、背景场核查、稳定化审计、原生P2、common-curved geometry、目标BC及两轮线性代数开发。
- 共享生产微泡模块的144个source/resource文件与189个active test文件，共333项，逐文件校验一致。没有整库stage共享工作树或恢复已删除的历史测试。
- 四轮原始监督job.log和EXECUTION.json，以及最新版已有完整evidence。服务器补充139个文件，共17,423,564 bytes，均与服务器SHA256一致。
- 必需的5.92 MB frozen curved-elevation NPZ、原P1 step10 checkpoint和对应VTU，以及本轮使用的BC/XML/PETSc选项文件。全量首系统矩阵转储未复制。
- 当前正式旧1565轨迹账本/summary/身份及既有必要可视化仍由原仓库文件保留。相对兼容链接按链接同步，没有复制成重复大文件。

包含清单共2,816项（不含本同步元数据与仓库继承文件），全都保持来源字节或symlink目标。远端已有、但不在独立算例index中的50项保持原样，见[继承文件](INHERITED_REMOTE_ONLY.json)；其它项目内容继承原Git树，不重新复制或更新。

## 大文件选择

新增/修改文件总量约143 MB，最大新增文件是21.28 MB的 `STABILIZATION_TERMS.npz`。这不是新压缩上传：直接使用已存在的无损科学NPZ。没有新增≥25 MB文件，没有启用LFS、切片、改写历史或force push。

**不上传86,869,515 bytes的重复 `STABILIZATION_TERMS.csv`。** 对全部271,473行及所有标识/数值列完成逐行比较，与已保留NPZ一致（包括NaN对应）。保留原导出代码 `terms.py` 和[验证记录](CSV_NPZ_EQUIVALENCE.json)，CSV原文件留在本地及服务器。原历史报告中的CSV链接/历史OUTPUT_SHA256仍指向这份刻意不发布的原文件；本轮不改写这些科学证据。NPZ和JSON结果可直接使用。

其余排除：虚拟环境、缓存、编译对象/二进制、完整第三方checkout/build、原PETSc矩阵大转储、冗余瞬态输出，以及旧清单中的大场数组、全壁表和旧动画副本。见[本轮排除](EXCLUDED_FILES.json)、[服务器大runtime目录项](SERVER_EVIDENCE.json)以及原 `EXCLUDED_LARGE_ARTIFACTS.json` / `GITHUB_SYNC_EXCLUDED.json`。原文件没有被删除，仍应按原路径与SHA取回；不会为取回数据自动重跑旧计算。

## 科学状态不变

P1同检查点解等价回归通过。P2 CPU Schur setup在任务RSS达到25.22 GiB时触发预声明24 GiB上限，未发生主机OOM或CUDA分配失败。停前保守余量9.90 GiB高于所需8.41 GiB，但内存仍增长、setup未完成；不能称完成后的安全门通过。

最终 `P2_PIPE_LINEAR_ALGEBRA_FEASIBLE=NO` / `STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT` 只限定本轮保护范围。完整P2 CFD未运行，曲面压力回填gate仍FAIL；不能由P1回归、工程检查或上传成功推断场精度通过。正式旧1565条仍含141截止、2失败，不是新固定灌注协议结果。

## 发布核验与复现边界

[SYNC_VALIDATION.json](SYNC_VALIDATION.json)验证包含文件SHA、当前报告输出manifest、改动范围、Python语法、JSON格式、新文件大小和凭据规则。保持根 `.gitattributes` 字节不变，不改科学文件行尾。已有249项工程检查及其失败/停止科学门均原样同步；同步阶段不重新求解。

配置保留原绝对路径和原环境契约，不宣称这是空白机器开箱即跑的安装包；源码和必要审计数据可核查，编译环境/大型原矩阵需按manifest恢复。GitHub不直接执行HTML，可下载相应目录后打开或直接查看PNG。

实际发布的commit SHA以本分支Git HEAD为准；push后在本地操作回执记录远端SHA校验，避免在提交内容中嵌入自身SHA。使用普通push创建新分支，并重新获取远端ref验证，结束后不启动任何科学任务。
