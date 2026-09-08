# 来源、复用与许可证

本地实际 Mirheo 库和源码为权威，未从同步分支覆盖本地文件。
源码提交、共享库哈希、精度、补丁和逐项函数/行号证据写入准备包及结果包的
`native_sdpd_contract.json` 与 `provenance.json`。

| 能力 | 复用来源 | 本次处理 |
| --- | --- | --- |
| 冻结物理目标与精确单位 | 第二阶段指定结果包的 `physical_case.json`、`unit_mapping.json` | 验证包清单哈希并原值读取；不按舍入值重建单位 |
| SI 换算、路径、哈希、JSON | `fluid_physics/units.py`、`common.py`、`vessel_geometry/io.py` | 直接调用 |
| GPU 超时、缓存与进程 | `fluid_physics/runner.py` | 增加共享预算成员注册读取；复用原 GPU 锁、monotonic 计费与自身进程组清理 |
| 统计、黏度、压力 | `fluid_physics/analysis.py`、`reporting.py` | 复用并修正候选分组、CSV 时间量化、温度稳定性、实际密度与压力覆盖边界 |
| 周期 GPU 采样循环 | `fluid_physics/gpu_worker.py` | 新比较 worker 改编相同两 rank/chunk 控制，增加原生 Density/SDPD、同步计时、原生密度快照和温度细分箱；旧 worker 原样保留 |
| 原生 SDPD 接线 | Mirheo `tests/sdpd/rest.py`、`double_poiseuille.py`、`rigid.py` | 使用同核同 rc 的 Density+SDPD，保留接口来源；未使用示例黏度替代目标物性 |
| 核函数及压力定义 | 本地 `sdpd.h`、`density*.h`、`pressure_EOS.h`、`stress_wrapper.h`、plugins 与 bindings | CPU 源码契约及少量真实快照核查；不改 C++/CUDA |
| 离线页面与 Windows 打开 | `test_code/review_fluid_physics.py`、`review_vessel_geometry.py` | 同 Plotly 基础与原打开函数，新增对照表和真实曲线 |
| 浏览器实测 | `test_code/check_fluid_physics_browser.cjs` | 复用 Windows Chrome/CDP 进程与截图方法，检查新页面交互 |

Mirheo 上游 MIT 许可证在每个 GPU attempt 中原样保存为 `MIRHEO_LICENSE.txt`。
上游源码 Copyright 注释不修改。CPU 核函数重算注明来源与数学公式，只用于观测核查。
两个旧用户工程没有据以覆盖整个工程的统一许可证，本次不附加这样的声明。
