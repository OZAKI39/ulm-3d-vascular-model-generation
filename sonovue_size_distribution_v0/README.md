# SonoVue size distribution V0

本项目只负责冻结 histogram 的连续 CDF 与逆 CDF 采样。

先读 [中文报告](SONOVUE_SAMPLER_REPORT.md)、[冻结合同](contracts/SONOVUE_SAMPLER_CONTRACT_V0.json) 和 [验证结果](validation/SONOVUE_SAMPLER_VALIDATION.json)。

已验证解释器：`/usr/bin/python3`；NumPy 1.26.4，Matplotlib 3.6.3（只用于 QA）。采样核心只依赖 NumPy 与 Python 标准库，不需要 SciPy。

```python
from src.sonovue_sampler import sample_diameters, SonoVueDistribution
diameters_um = sample_diameters(n=1000, seed=20260915)  # demo only
distribution = SonoVueDistribution()
d50 = distribution.inverse_cdf(0.5)
```

```bash
/usr/bin/python3 -B -m unittest discover -s tests -v
sha256sum -c SHA256SUMS
```

原始冻结输入及合同不覆盖；CLI 拒绝覆盖已有 population。`src/validate_sampler.py` 保留本次真实验证流程，已封存的 validation 不会自动覆盖；若需再次审查，使用新的独立副本，保留本次证据。`src/finalize_contract.py` 只封存已有验证，不抽样、不启动后续任务。

所有现有 population 都是 DEMO / VALIDATION，正式数量与浓度尚未定义。NEXT_STAGE=MINIMAL_LAMMPS_PARTICLE_ENGINE，仅为下一阶段名称，未启动。
