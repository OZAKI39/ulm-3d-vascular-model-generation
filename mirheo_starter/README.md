# Mirheo starter：代码与实测数据归档

本目录同步自 `/home/lzy/projects/mirheo_starter`，分支为
`sync/mirheo-starter-20260908`。2026-09-09 本次同步包含 **18,425 个源文件、
479,123,181 字节**，保留代码、配置、数据、原始运行、账本、报告与截图的原始字节和权限。
文件级 SHA256 与排除清单见 [SYNC_MANIFEST.json](SYNC_MANIFEST.json)，
同步验证见 [SYNC_VALIDATION.json](SYNC_VALIDATION.json)。

## 最新交付

固定参数无驱动 SDPD 实验已完成 **400000 步 / 12.662309354 µs**。
正式区间为 `t*∈(0.30,0.40]`，包含 500 个样本。
科学状态为 **`STATIONARITY_OR_SAMPLING_INCONCLUSIVE`**：温度两半均值差约
0.790%，超过事前 0.5% 门槛；压力和最近邻结构也未通过筛选。
尚未建立有效稳态块数、稳定温度均值或 CI，`selection=null`。
程序完整运行不等于液体全部物性通过。

- [实测中文报告](test_code/outputs/sdpd_equilibration/result_26300860927d8b87/fcb8821bfdea/report_zh.md)
- [离线可视化 HTML](test_code/outputs/sdpd_equilibration/result_26300860927d8b87/fcb8821bfdea/sdpd_equilibration_review.html)：下载完整目录后直接用浏览器打开；无需 GPU。
- [真实分析数据包](data/sdpd_equilibration/result_26300860927d8b87/)
- [完整原始轨迹、快照、日志与执行记录](runs/sdpd_equilibration/unforced_plateau_20260909/longer_unforced_thermal_plateau/)
- [原交付记录与验证哈希](test_code/outputs/sdpd_equilibration/execution_20260909T114110_487235Z/delivery_record.json)
- [实验入口及细节](test_code/README_sdpd_equilibration.md)

原实验的 200 项 CPU 测试、15 项真实浏览器检查均通过，人工验收仍为 PENDING。
本次 GPU 任务实收 528.365623583 s；用户追加 493 s 后总授权为 4093 s，
累计用量 4021.284952471 s，剩余 71.715047529 s。
授权记录与所有历史用量一同归档；仓库复制不会产生新预算，也不允许重复记入追加授权。

首次 JSON 导出失败后的原始证据保留；仅在 CPU 端修复 NumPy 布尔值序列化并重新分析，
没有重跑 GPU。`data/sdpd_equilibration/result_aa25ea397a1f12de/` 为不完整失败导出，
不能作为有效结果；应使用以上明确交付的数据包。历史输出中的旧状态不覆盖最新交付。

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
200 项测试仍需要原环境、固定的外部来源，以及本地 Mirheo 源码和编译库。
归档数据不是已经完成环境迁移的可执行 GPU 安装包。

CPU 依赖见 [requirements-cpu.txt](requirements-cpu.txt)。在具备原环境的 WSL 中：

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m unittest discover -s test_code -p 'test_*.py' -v
.venv/bin/python -B -m test_code.review_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration.yaml --open
```

`--execute` 仍受原预算、精确计划与已有任务缓存约束。本次同步未新增 GPU 实验，
未修改液体物性结论，也未将人工验收状态改为通过。
