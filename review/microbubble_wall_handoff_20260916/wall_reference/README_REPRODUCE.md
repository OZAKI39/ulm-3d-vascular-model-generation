# Wall reference audit 复现说明

本归档状态为 **PASS_WITH_LIMITATIONS**。先阅读 `WALL_HYDRODYNAMICS_REFERENCE_AUDIT_REPORT.md` 和 `PRIMARY_REFERENCE_SELECTION.json`；不要把 finite 或官方测试通过理解为完整矩阵正确。

## 已完成的执行

全部计算在本地 WSL CPU 执行，线程限制为 `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`，并设置 `CUDA_VISIBLE_DEVICES=''`。没有访问 Vast，没有修改现有科研环境或生产代码。

工作目录：

```text
/home/lzy/projects/wall_hydrodynamics_reference_audit_20260916_115322
```

结果目录：

```text
/home/lzy/projects/compre_output/wall_hydrodynamics_reference_audit/20260916_115322
```

`upstream/RigidMultiblobsWall` 和 `upstream/pystokes` 是冻结官方 checkout；实际构建在 `build/`，环境在 `envs/wallref_rmbw`、`envs/wallref_pystokes`。环境初始和最终 pip manifest 同时保留。RMBW conda 的精确构建规格位于 `provenance/wallref_rmbw_conda_explicit.txt`；Python 3.12 环境最终规格见 `provenance/wallref_pystokes_pip_freeze_FINAL.txt`。

## 校验和审阅

在结果目录执行：

```bash
sha256sum --check --quiet SHA256SUMS
cat FINAL_TERMINAL_SUMMARY.txt
```

`SHA256SUMS` 包含该目录内除其自身以外的全部最终文件；`CODE_SHA256SUMS` 单独列出 audit Python 源码。最终校验回执保存在工作目录的 `RESULT_SHA256_VERIFICATION.json`，避免 manifest 对自身/回执的循环引用。上游源码 tar 是包含源码与许可证的参考工具归档，没有把 GPL 源码复制成生产插件。

结果链条为：

```text
冻结 source + SonoVue contract + WALL_REFERENCE_CONVENTION.json
  -> raw/{RMBW_SPHERE,RMBW_LUBRICATION,PYSTOKES}.jsonl
  -> audit_wall_reference.py
  -> CSV / HDF5 / matrix audits / classical comparisons
  -> visualize_wall_audit.py
  -> 8 QA PNG + VISUALIZATION_PROVENANCE.json
  -> finalize_audit_evidence.py
  -> write_final_report.py
  -> report / decision / final summary / SHA256SUMS
```

378 个主 sweep 矩阵和 16 个历史表格诊断都是首次输出。`FIRST_OUTPUT_IDENTITY.json` 记录其最终核查前的字节身份；collector 拒绝覆盖已有文件。旧 sphere 首次导入的两条警告和所有科学失败均保留。finalizer 最初错误地要求所有 import 无警告，已改为保留并报告警告；未因该审计脚本修正重跑 reference。

## Reference API

两个适配器均返回 `(matrix_SI, metadata)`，行顺序 Vx,Vy,Vz,Ωx,Ωy,Ωz，列顺序 Fx,Fy,Fz,Tx,Ty,Tz。

```python
# 必须在 RMBW 独立 Python 3.8 环境内导入。
from rmbw_reference_adapter import get_wall_mobility
M, metadata = get_wall_mobility(radius_m, gap_m, viscosity_pa_s,
                                implementation="lubrication")
```

为保留原正式 collector 行为，RMBW adapter 原默认参数仍为 `sphere_semianalytical`；**使用选定主参考时必须显式指定 `implementation="lubrication"`**。collector 对两个入口均显式指定该参数。不要使用旧默认入口的完整矩阵作为生产基准。

```python
# 在独立 Python 3.12 环境内导入；完整矩阵存在已知互易性问题。
from pystokes_reference_adapter import get_wall_mobility
M, metadata = get_wall_mobility(radius_m, gap_m, viscosity_pa_s)
```

RMBW 原始 Lubrication API 返回 excess R。wrapper 在固定 native 单位下加回 isolated bulk、求六个单位负载响应、转换 SI。每个 case 的 native R/M、单位因子、cutoff 状态均写入 raw metadata，独立 finalizer 再次核对转换。没有加入模型修复或符号补丁。

## 在当前环境独立复现数值

以下命令**仅供后续显式复现**，本次交付没有额外执行。新建目录，避免覆盖当前证据；无需更改冻结 contract。

```bash
WALLREF_ORIGINAL=/home/lzy/projects/compre_output/wall_hydrodynamics_reference_audit/20260916_115322
WALLREF_WORK=/home/lzy/projects/wall_hydrodynamics_reference_audit_20260916_115322
WALLREF_REPLAY=$(mktemp -d /home/lzy/projects/wallref_replay.XXXXXXXX)
mkdir -p "$WALLREF_REPLAY"/raw "$WALLREF_REPLAY"/provenance "$WALLREF_REPLAY"/validation "$WALLREF_REPLAY"/visualization
cp -a "$WALLREF_ORIGINAL/src" "$WALLREF_REPLAY/src"
cp "$WALLREF_ORIGINAL/WALL_REFERENCE_CONVENTION.json" "$WALLREF_REPLAY/"
cp "$WALLREF_ORIGINAL/provenance/TASK_PATHS.json" "$WALLREF_REPLAY/provenance/"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=''
"$WALLREF_WORK/envs/wallref_rmbw/bin/python" -B "$WALLREF_REPLAY/src/collect_package_results.py" RMBW_SPHERE
"$WALLREF_WORK/envs/wallref_rmbw/bin/python" -B "$WALLREF_REPLAY/src/collect_package_results.py" RMBW_LUBRICATION
"$WALLREF_WORK/envs/wallref_pystokes/bin/python" -B "$WALLREF_REPLAY/src/collect_package_results.py" PYSTOKES
"$WALLREF_WORK/envs/wallref_rmbw/bin/python" -B "$WALLREF_REPLAY/src/compare_oneill_table.py" RMBW_LUBRICATION
"$WALLREF_WORK/envs/wallref_pystokes/bin/python" -B "$WALLREF_REPLAY/src/compare_oneill_table.py" PYSTOKES
"$WALLREF_WORK/envs/wallref_pystokes/bin/python" -B "$WALLREF_REPLAY/src/audit_wall_reference.py"
"$WALLREF_WORK/envs/wallref_pystokes/bin/python" -B "$WALLREF_REPLAY/src/visualize_wall_audit.py"
```

模型结果应按矩阵数值、状态及已知异常比较；运行 wall_seconds、日志时间和环境警告不是 byte-identical 的物理结果要求。不要用再次执行来替换第一次输出。

只重做后处理时，将 `raw/` 复制到新的 replay 目录并跳过五条 collector/history 命令。完整源码/baseline finalizer 还需要 `provenance/`、package provenance 和可访问的原始 baseline 路径；它不是无原始 WSL 依赖的远程运行程序。

## 在另一台机器重建

归档 tar 的顶层直接是各仓库文件，需要分别解压到独立空目录。tar 不包含 `.git`，commit/gitlink 来源见 INITIAL_SOURCE_IDENTITY；源码内容可逐项用对应 SOURCE_SHA256SUMS 核对。若需要 git 终审，则从官方仓库取得并 checkout 指定 commit。

用 RMBW 的 conda explicit 文件创建独立 Python 3.8 环境。编译命令的原始 argv 在 `provenance/RMBW_LUB_BUILD_COMMAND.json`，链接依赖见 `RMBW_LUB_LINKAGE.txt`。在复制件中替换工作目录和环境前缀；`__FILENAME__` 必须指向该机器冻结源码的 `Lubrication` 目录，以读取原表格。实际编译使用系统 Eigen 3.4 头文件，而不是 conda 内另装的 Eigen 5。不得编辑上游 `.cc` 或公式。

PyStokes 使用独立 Python 3.12 venv，按最终 pip manifest 中的精确版本安装依赖，再从冻结 source 副本以 `pip install --no-build-isolation` 构建；manifest 中本机 `file://` 来源需指向该副本。官方短测试可用该环境的 `python -B upstream/pystokes/tests/test_short.py` 运行，13 个测试结果见 `logs/PYSTOKES_OFFICIAL_SHORT_TEST.log`。这些步骤用于重现已冻结 2.3.3，不追随仓库新版本。

在新的结果副本 `provenance/TASK_PATHS.json` 中更新 `work`；wrapper 从该键读取 source/build 路径。这是复现路径配置，不是物理 contract 变更。本归档未声称可跨 Python/编译器平台得到逐字节相同的浮点结果。

## 停止状态

`HUMAN_VISUAL_REVIEW=PENDING`。不自动修正 PyStokes、旧 RMBW sphere，或编写任何正式 wall module。下一步为用户审阅，微泡 wall hydrodynamics / exclusion / adhesion / local-plane / curved STL 均保持 PENDING，RBC=OFF。
