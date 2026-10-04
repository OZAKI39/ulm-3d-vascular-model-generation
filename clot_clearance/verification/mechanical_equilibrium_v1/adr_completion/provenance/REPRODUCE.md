# 复现和文件身份

工作树：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance`
分支：`clot-clearance/energy-regularization-20260930`
基线提交：`f0e891686b9b107510ab3292aa27ce83638b3339`
未自动提交或推送本轮新增文件。

Python：`/home/lzy/projects/computation_examples/_shared/python_environment/bin/python`
依赖版本与硬件见 `../HARDWARE.json`，所有数值计算沿用现有环境。
PDF 核对用 pypdf 及只安装于忽略目录 `build_mechanical_equilibrium_v1/adr_completion/pdf_tools` 的 PyMuPDF。

环境变量：

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMBA_CACHE_DIR=build_mechanical_equilibrium_v1/adr_completion/reproduction_numba
export MPLCONFIGDIR=build_mechanical_equilibrium_v1/adr_completion/reproduction_matplotlib
```

在**新工作副本**中依次运行（脚本拒绝覆盖已有交付目录；不要删除此交付来重跑）：

```text
python -B scripts/verify_adr_completion.py small
python -B scripts/verify_adr_completion.py erosion
python -B -m pytest tests -q -p no:cacheprovider --basetemp=<new test directory>
python -B scripts/audit_adr_completion.py
python -B scripts/review_adr_completion.py
```

首次无保护分支 probe 的程序、输入和原诊断均在 `../initial_small_probe/`；后来用原代码精确复现一次以补存完整停止状态，
其 SUMMARY 逐字段完全相同。该失败不是正式 guarded campaign，也没有发生损伤更新。
正式 source 使用 guards 的依据详见 `../DAMPING_GUARD_ADDENDUM.md`。

所有旧78项与前置15项测试未修改。新19项测试的输出保存在 `../test_artifacts/`。
每个正式job有单独log、输入config、源码快照和SHA；结果保存全部平衡阶段与全部35次damage jump。
每次平衡归零的是虚拟速度，既不是治疗时间速度，也不是脱落微粒运动轨迹。
末次损伤跳变后再平衡的状态用于最终图，不再推进曝光或损伤。

保存的 NPZ/CSV 缺失数值以NaN表示“该列此步不适用”，并非科学计算的非有限结果。
JSON中的 physical velocity/time/kinetic energy 为 null，均没有伪装为已计算的物理量。
所有图直接使用保存数据，不平滑，不改CSV，不放大位移。

范围限制：完整锚定拓扑、静态mean flow、当前未标定材料及损伤演示参数。
没有将旧1565条微泡轨迹或 BraVa 正式CFD重标记为本轮结果。
