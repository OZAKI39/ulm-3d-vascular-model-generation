# Mirheo starter：代码与实测数据归档

本目录同步自 `/home/lzy/projects/mirheo_starter`，分支为
`sync/mirheo-starter-20260908`。2026-09-09 本次同步包含 **19,022 个源文件、
546,902,551 字节**，保留代码、配置、数据、原始运行、账本、报告与截图的原始字节和权限。
文件级 SHA256 与排除清单见 [SYNC_MANIFEST.json](SYNC_MANIFEST.json)，
同步验证见 [SYNC_VALIDATION.json](SYNC_VALIDATION.json)。

## 最新交付

最新交付为 **`RESTART_NOT_VALIDATED / RESTART_DIAGNOSTIC_FAILED`**。
已实际执行三进程非零 dt 恢复对照：A 连续 4000 步完成；B 推进 2001 步，
归档第 2000 步的真实状态；B 恢复进程在原生读取时退出，未得到恢复后推进前状态。
直接原因是 `saved_forces` 被写成 `Other / Force`、HDF 形状 `[4096,1]`，
原生读取器不支持该格式。`saved_stresses` 的 Tensor6 格式正常。
归档位置、速度和粒子 ID 与 B 自身保存边界完全一致，但不能据此宣布完整恢复通过。
原生 SDPD 相互作用 RNG 未持久化的问题也仍存在；未启动后续长实验，未自动重试。

本次同步包含新的 Python 通道格式预检、固定执行计划重分析、完整失败证据、
按 ID 的保存状态比较、实际状态与预算修复、独立短诊断曲线及保护测试。
物理参数、单位、原生库和冻结 GPU 执行脚本保持不变。

- [最新中文报告](data/sdpd_equilibration_extended/result_2f3922c3ebd6ca0d/report_zh.md)
- [最新离线可视化 HTML](test_code/outputs/sdpd_equilibration_extended/result_2f3922c3ebd6ca0d/9c15803c55e3/sdpd_equilibration_extended_review.html)：下载完整目录后直接打开；无需 GPU。
- [最新分析数据包](data/sdpd_equilibration_extended/result_2f3922c3ebd6ca0d/)
- [三进程原始数据、冻结代码、checkpoint 与账本](runs/sdpd_equilibration_extended/fixed_late_restart_20260909/)
- [最新交付记录及验证哈希](test_code/outputs/sdpd_equilibration_extended/execution_20260909T153829_890721Z/delivery_record.json)
- [259 项 CPU 测试记录](test_code/outputs/sdpd_equilibration_extended/execution_20260909T153829_890721Z/validation_after_fixes/CPU_validation.json)
- [20 项真实浏览器检查](test_code/outputs/sdpd_equilibration_extended/result_2f3922c3ebd6ca0d/9c15803c55e3/browser_final/browser_checks.json)

液体长观察仍止于 **400000 步 / t*=0.40 / 12.662309354 µs**。
旧正式区间 `(0.30,0.40]` 有 500 个样本；温度两半差约 0.790%，超过原 0.5% 门槛，
补充分析仍支持冷却。压力瞬时波动强、持续总压力漂移未证实；核密度及最近邻分布仍变化。
旧科学结论 `STATIONARITY_OR_SAMPLING_INCONCLUSIVE` 保留；新正式窗口 `(0.60,0.80]`
没有实测样本。本次短恢复诊断不作为液体平衡证据，`selection=null`，人工验收 PENDING。
[旧长实验报告](test_code/outputs/sdpd_equilibration/result_26300860927d8b87/fcb8821bfdea/report_zh.md)
和[完整旧轨迹](runs/sdpd_equilibration/unforced_plateau_20260909/longer_unforced_thermal_plateau/)均保留。

本次恢复诊断实收 **11.711100827 s**。追加 1379 s 的用户授权已登记，
总授权为 **5472 s**，累计用量 **4032.996053298 s**，范围内剩余 **1439.003946702 s**。
余额中原基额剩余 71.715047529 s，新用途扩展剩余 1367.288899173 s；旧 493 s 扩展已耗尽。
1360 s 条件长段预算未支出，预算余额不会解除恢复限制。
授权记录与所有历史用量一同归档；仓库复制不会产生新预算，也不允许重复记入追加授权。

历史失败证据均保留。`data/sdpd_equilibration/result_aa25ea397a1f12de/` 是旧的不完整导出；
`data/sdpd_equilibration_extended/preparation_779eb32f8aae47f0/` 的执行状态和重复预算请求
存在已修复的 Python 汇报错误。它们不覆盖上方明确交付的最新结果。

## 同步范围

| 目录 | 内容 |
| --- | --- |
| `py_scripts/`、`scripts/`、`test_code/` | 现有源码、配置、自动测试、复核入口 |
| `data/` | 已验收几何、单位与工况、液体对照、SDPD 诊断及稳定性数据包 |
| `runs/` | 原始轨迹、同相位快照、冻结执行代码、日志、完整预算与授权账本 |
| `test_code/outputs/` | 交付记录、历史及最终报告、离线 HTML、截图、CPU 与浏览器验证记录 |
| `metadata/`、`logs/` | 原 WSL 环境、原生版本、兼容补丁及安装/运行日志 |

本次包含历史与失败证据，没有清理或改写它们。虚拟环境、第三方 Mirheo 检出/构建、
下载目录、字节码、运行锁和临时浏览器配置不纳入 Git；精确排除项见同步清单。
所有纳入的源文件均小于 10 MB，直接保存为 Git 对象；不使用 Git LFS 指针。
`.gitattributes` 禁用换行转换，防止克隆后的字节与原 SHA256 不符。
上游 Mirheo 许可证副本见 [MIRHEO_LICENSE.txt](MIRHEO_LICENSE.txt)。

## 代码入口

| 内容 | 入口与说明 |
| --- | --- |
| 已验收血管几何与标签迁移 | `py_scripts/import_vessel_geometry.py`；[说明](test_code/README_vessel_geometry.md) |
| 工况、单位与 DPD 标定 | `py_scripts/prepare_fluid_physics.py`、`py_scripts/calibrate_dpd_fluid.py`；[说明](test_code/README_fluid_physics.md) |
| DPD/SDPD 纯液体对照 | `py_scripts/compare_dpd_sdpd.py`；[说明](test_code/README_fluid_model_comparison.md) |
| SDPD 原因排查 | `py_scripts/diagnose_sdpd.py`；[说明](test_code/README_sdpd_diagnostics.md) |
| 固定参数持续无驱动演化 | `py_scripts/run_sdpd_equilibration.py`；[说明](test_code/README_sdpd_equilibration.md) |
| checkpoint 恢复诊断与后段计划 | 同一入口配合 `py_scripts/sdpd_equilibration_extended.yaml`；[最新报告](data/sdpd_equilibration_extended/result_2f3922c3ebd6ca0d/report_zh.md) |
| 离线复核 | `test_code/review_*.py`、`test_code/check_*_browser.cjs` |

## 核对与运行位置

在本目录运行以下 CPU 命令，可以逐文件核对同步清单，不启动模拟：

```bash
python3 - <<'PYTHON'
import hashlib, json
from pathlib import Path
manifest = json.loads(Path('SYNC_MANIFEST.json').read_text())
for name, expected in manifest['files'].items():
    path = Path(name)
    assert path.stat().st_size == expected['size'], name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected['sha256'], name
print('PASS:', len(manifest['files']), 'files')
PYTHON
```

源码、YAML、HTML 和 provenance 中原 WSL 路径及哈希保持原样，不能批量替换路径后
仍声称匹配原始证据。浏览已归档 HTML 不需要 Mirheo；重新运行完整诊断、GPU 或全部
259 项测试仍需要原环境、固定的外部来源，以及本地 Mirheo 源码和编译库。
归档数据不是已经完成环境迁移的可执行 GPU 安装包。

CPU 依赖见 [requirements-cpu.txt](requirements-cpu.txt)。在具备原环境的 WSL 中：

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m unittest discover -s test_code -p 'test_*.py' -v
.venv/bin/python -B -m test_code.review_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration_extended.yaml --open
```

`--execute` 仍受原预算、精确计划与已有任务缓存约束。本次同步未新增 GPU 实验，
未修改液体物性结论，也未将人工验收状态改为通过。
