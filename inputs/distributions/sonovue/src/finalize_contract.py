#!/usr/bin/env python3
"""Freeze the completed sampler audit; no further populations are generated."""
from pathlib import Path
import ast
import hashlib
import json
import platform

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def main():
    contract_path = ROOT / 'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json'
    if contract_path.exists():
        raise FileExistsError('Sampler contract already frozen; do not overwrite')
    v = json.loads((ROOT / 'validation/SONOVUE_SAMPLER_VALIDATION.json').read_text())
    assert v['status'] == 'PASS' and all(s == 'PASS' for s in v['checks'].values())
    assert v['sampler_source_sha256'] == sha(ROOT / 'src/sonovue_sampler.py')
    assert v['validation_source_sha256'] == sha(ROOT / 'src/validate_sampler.py')
    assert sha(ROOT / 'input/FROZEN_SONOVUE_HISTOGRAM.csv') == v['histogram_sha256']
    assert sha(Path('/home/lzy/projects/FROZEN_SONOVUE_HISTOGRAM.csv')) == v['histogram_sha256']
    source = (ROOT / 'src/sonovue_sampler.py').read_text()
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in ('sample_diameters', 'sample_with_provenance', 'inverse_cdf', '_inverse_with_indices'):
            assert not any(isinstance(x, (ast.For, ast.While, ast.ListComp, ast.GeneratorExp)) for x in ast.walk(node)), node.name
    assert 'np.random.seed(' not in source and 'scipy' not in source
    write_json(ROOT / 'provenance/SOURCE_STATIC_AUDIT.json', {
        'status': 'PASS', 'vectorized_sampling_path': True, 'per_particle_Python_loops': False,
        'global_numpy_seed_used': False, 'scipy_or_parametric_fit': False,
        'sampler_source_sha256': v['sampler_source_sha256'],
        'inspection_scope': 'AST of sample/inverse methods; CSV serialization and bin-table parsing are outside the sampling timer',
    })

    contract = {
        'contract_name': 'SONOVUE_SAMPLER_CONTRACT_V0', 'status': 'FROZEN_VALIDATED_PASS',
        'primary_source': 'Kotopoulis et al., Pharmaceutics 2022 Figure 3A',
        'primary_source_doi': '10.3390/pharmaceutics14010098',
        'primary_source_url': 'https://www.mdpi.com/1999-4923/14/1/98',
        'population': 'fresh SonoVue', 'weighting': 'count-normalized / number-weighted',
        'histogram_file': 'input/FROZEN_SONOVUE_HISTOGRAM.csv', 'path_base': 'project root',
        'histogram_sha256': v['histogram_sha256'], 'histogram_original_modified': False,
        'bin_count': v['bin_count'], 'zero_probability_bins': v['zero_probability_bins'],
        'bin_width_um': v['bin_width_um'], 'bin_width_min_um': v['bin_width_min_um'], 'bin_width_max_um': v['bin_width_max_um'],
        'original_bin_width_definition_um': 0.1, 'bin_width_absolute_check_tolerance_um': 1e-10,
        'support_min_um': v['support_min_um'], 'support_max_um': v['support_max_um'],
        'support_authority': 'First diameter_low_um and last diameter_high_um in the frozen CSV; never a separate hard-coded support',
        'representation': 'continuous empirical CDF', 'cdf_type': 'piecewise linear', 'pdf_type': 'piecewise constant',
        'representation_identifier': 'CONTINUOUS_PIECEWISE_LINEAR_EMPIRICAL_CDF',
        'sampling_method': 'inverse-CDF', 'within_bin_assumption': 'uniform',
        'explanation': 'Interpret each original 0.1 um histogram bin as a constant PDF with its empirical probability mass; its CDF is linear. This is not a new parametric fit and uses no smoothing or assumed distribution family.',
        'benefits': ['continuous diameter output', 'preserve empirical bin masses', 'no extra smoothing', 'no assumed distribution family'],
        'normalization_method': 'p_i divided by the common math.fsum of input probability; no per-bin reweighting',
        'INPUT_PROBABILITY_SUM': v['INPUT_PROBABILITY_SUM'], 'NORMALIZED_PROBABILITY_SUM': v['NORMALIZED_PROBABILITY_SUM'],
        'cdf_float64_endpoint_policy': 'Accumulate normalized probabilities in source order; set last positive cumulative endpoint and any trailing zero endpoints to exactly 1. Zero masses remain exactly zero.',
        'cdf_mass_vs_normalized_probability_max_abs_roundoff': v['cdf_mass_vs_normalized_probability_max_abs_roundoff'],
        'inverse_interval_denominator': 'cdf_high - cdf_low; agrees with normalized probability to recorded float64 roundoff',
        'inverse_plateau_convention': 'Right-side: first positive interval with cdf_high > u for 0<=u<1. Exact plateau probability maps to next positive low edge; inverse(0)=first positive low; inverse(1)=last positive high.',
        'sample_interval_convention': 'Positive [low,high) intervals; exact floating-point high-edge rounding for u<1 is guarded with nextafter(high,low). Shared boundaries have zero continuous mass.',
        'zero_probability_bin_policy': 'PRESERVED; CDF flat; never sampled as an interval',
        'zero_probability_bin_count': v['zero_probability_bins'],
        'gap_policy': 'Any non-overlapping gap remains zero density; no adjacent bin widths changed',
        'parametric_fit': 'NO', 'Gaussian_fit': 'NO', 'lognormal_fit': 'NO', 'KDE': 'NO', 'smoothing': 'NO',
        'bin_merging': 'NO', 'tail_truncation': 'NO', 'tail_extension': 'NO', 'volume_weighted': 'OUT_OF_SCOPE',
        'reference_mean_diameter_um': 2.51, 'reference_mean_role': 'INFORMATIONAL_ONLY', 'reweight_to_2p51': 'NO',
        'target_continuous_mean_um': v['TARGET_MEAN_DIAMETER_UM'], 'target_quantiles_um': v['target_quantiles_um'],
        'random_generator': 'numpy.random.default_rng', 'observed_bit_generator': v['environment']['bit_generator'],
        'seed_policy': 'frozen_per_case', 'global_distribution_seed': None,
        'population_required_metadata': ['seed', 'N', 'histogram_sha256', 'sampler_version'],
        'population_metadata_format': 'CSV plus mandatory adjacent .metadata.json bound to CSV SHA256',
        'sampling_api': 'sample_diameters(n, seed, histogram=DEFAULT_HISTOGRAM)',
        'sampler_version': v['sampler_version'], 'sampler_source_file': 'src/sonovue_sampler.py',
        'sampler_source_sha256': v['sampler_source_sha256'], 'dtype': 'float64', 'rng_uniform_domain': '[0,1)',
        'reproducibility': {'status': 'PASS', 'same_histogram_n_seed': 'byte-identical arrays and CSV in two independent processes',
                            'numpy_version_tested': v['environment']['numpy'], 'numpy_major_tested': 1,
                            'cross_version_or_cross_architecture_guarantee': 'NOT_CLAIMED; record exact environment and output hash per case'},
        'validation_plan_file': 'contracts/VALIDATION_PLAN_V0.json', 'validation_plan_sha256': v['validation_plan_sha256'],
        'validation_file': 'validation/SONOVUE_SAMPLER_VALIDATION.json', 'validation_sha256': sha(ROOT / 'validation/SONOVUE_SAMPLER_VALIDATION.json'),
        'cdf_file': 'validation/SONOVUE_EMPIRICAL_CDF.csv', 'cdf_file_sha256': sha(ROOT / 'validation/SONOVUE_EMPIRICAL_CDF.csv'),
        'cdf_table_layout': 'Two rows per original bin (low/high); repeated shared edges are intentional; all zero bins remain present',
        'formal_simulation_population_created': False, 'formal_population_size': 'PENDING_FUTURE_SIMULATION_CASE_CONTRACT',
        'demo_population_role': 'SAMPLER_DEMO_ONLY', 'validation_population_role': 'SAMPLER_VALIDATION_ONLY',
        'performance_tuning': 'OUT_OF_SCOPE',
        'scientific_scope': 'Sampler correctness conditional on the frozen manually digitized histogram and uniform-within-bin assumption; no new validation of digitization or transport physics',
        'MICROBUBBLE_TRANSPORT_MODEL': 'PENDING', 'LAMMPS_PARTICLE_ENGINE': 'PENDING', 'PALABOS_LAMMPS_COUPLING': 'PENDING',
        'MICROBUBBLE_WALL_MODEL': 'PENDING', 'MICROBUBBLE_ADHESION': 'PENDING', 'RBC': 'OFF',
        'NEXT_STAGE': 'MINIMAL_LAMMPS_PARTICLE_ENGINE', 'NEXT_STAGE_STARTED': False,
    }
    write_json(contract_path, contract)
    q, sq = v['target_quantiles_um'], v['validation_quantiles_um']
    gate_rows = '\n'.join(f"| {name} | {status} |" for name,status in v['checks'].items())
    report = f'''# SonoVue 连续经验 CDF 与逆 CDF 采样器 V0

**SONOVUE_CONTINUOUS_SAMPLER_V0 = PASS**。本轮只验证尺寸分布采样，生成 1000 个 demo 样本与 100000 个验证样本，均不是正式 simulation population。下一阶段标识为 `MINIMAL_LAMMPS_PARTICLE_ENGINE`，本轮没有启动它。

## 输入与科学定义

来源标识为 [Kotopoulis et al., Pharmaceutics 2022, Figure 3A](https://www.mdpi.com/1999-4923/14/1/98)，DOI: 10.3390/pharmaceutics14010098。冻结对象为 fresh SonoVue、count-normalized / number-weighted 分布。文章引用仅用于标注来源；本轮没有重新数字化 Figure 3A，也没有导入另一套概率。

正式输入 `/home/lzy/projects/FROZEN_SONOVUE_HISTOGRAM.csv` 原样复制到 [input/](input/FROZEN_SONOVUE_HISTOGRAM.csv)，两者 SHA256 均为 `{v['histogram_sha256']}`。CSV 的中心、low/high 和 probability 是唯一输入；附带 cdf、digitized_frequency_percent 等列不参与计算。

共有 {v['bin_count']} 个原始 bin，含 {v['zero_probability_bins']} 个零概率 bin，均保留。support 从 CSV 读取为 **{v['support_min_um']}–{v['support_max_um']} µm**。浮点读取的 bin 宽度范围为 {v['bin_width_min_um']:.17g}–{v['bin_width_max_um']:.17g} µm，中位数 {v['bin_width_um']:.17g} µm；符合冻结的约 0.1 µm 定义，允许绝对浮点差 1e-10 µm。

输入概率和 `{v['INPUT_PROBABILITY_SUM']!r}`，只除以这个共同总和后为 `{v['NORMALIZED_PROBABILITY_SUM']!r}`。没有改变相对概率、合并 bin、平滑、Gaussian/lognormal/KDE/其他拟合、体积加权、截尾或延长尾部。

## 连续分布及边界约定

每个原始 bin 的 PDF 为 `p_i / (high_i - low_i)`，CDF 从 `cdf_low` 线性增至 `cdf_high`；support 两端严格为 0 和 1。它是在 uniform-within-bin 假设下解释现有 histogram，不是新的参数分布。CDF 的水平段保持原宽度。

`diameter = inverse_cdf(U)`，使用 `numpy.random.default_rng(seed)` 的一次向量化 U 抽样。以 `searchsorted(..., side="right")` 选择首个 `cdf_high > U` 的正概率区间，随后解析求逆；没有另抽第二个随机数。精确命中水平段概率时选择下一个正概率 bin 的 low；诊断 API 的 U=0/1 分别返回首/末正概率 bin 的外端点，正式抽样 U 属于 [0,1)。

浮点结果若在 U<1 时恰好舍入为 high，以 `nextafter(high, low)` 保持在正概率半开区间内。这只是一个 ULP 的端点表示保护；不改变科学 support 或尾部定义。原概率与累积端点差分的最大浮点差为 {v['cdf_mass_vs_normalized_probability_max_abs_roundoff']:.17g}；零 bin 的差分仍严格为 0。

[SONOVUE_EMPIRICAL_CDF.csv](validation/SONOVUE_EMPIRICAL_CDF.csv) 为每个 bin 保留 low/high 两行，共 {2*v['bin_count']} 行，含 diameter/cdf、原区间、概率质量、PDF、CDF 低/高值。相邻共享端点重复是可核查布局，不是重复概率质量。

## 真实验证结果

| 检查 | 状态 |
|---|---|
{gate_rows}

8 项解析/边界/非法输入单元测试通过，其中还检查无全局 RNG 状态修改、非必要列被忽略、原件 hash 与只读数组；见 [单元测试记录](provenance/UNIT_TEST_EXECUTION.json)。没有修改 frozen validation plan 的阈值。

逆一致性使用 {v['inverse_N']} 个 U、seed={v['inverse_seed']}。`max|cdf(inverse_cdf(U))-U| = {v['INVERSE_CDF_MAX_ERROR']:.17g}`；另用逐原始 bin 概率积分的独立向量公式，误差为 {v['independent_integral_inverse_max_error']:.17g}；所有 CDF 边界概率误差为 {v['inverse_boundary_max_error']}. 均低于 1e-10。

100000 个验证样本 seed=20260915，全部有限且在原 support 内；7 个零概率 bin 合计 **0 个样本**。逐原 bin 重建详见 [BIN_REPRODUCTION.csv](validation/BIN_REPRODUCTION.csv)。

| 指标 | 实测 | 冻结 smoke 门槛 |
|---|---:|---:|
| max absolute bin probability error | {v['MAX_BIN_PROBABILITY_ERROR']:.17g} | ≤0.01 |
| RMSE bin probability | {v['RMSE_BIN_PROBABILITY']:.17g} | 仅报告 |
| KS-like max CDF deviation | {v['MAX_CDF_DEVIATION']:.17g} | ≤0.02 |
| absolute mean difference / µm | {v['mean_absolute_difference_um']:.17g} | ≤0.05 |

KS-like 指标对全部排序样本的 ECDF 左右极限取最大差：`max(i/N-F(x_i), F(x_i)-(i-1)/N)`。这是采样实现 smoke，没有计算 p-value，也不作为出版级拟合检验。

| 诊断 / µm | Target continuous distribution | Validation sample |
|---|---:|---:|
| mean | {v['TARGET_MEAN_DIAMETER_UM']:.17g} | {v['VALIDATION_SAMPLE_MEAN_UM']:.17g} |
| d10 | {q['d10']:.17g} | {sq['d10']:.17g} |
| d50 | {q['d50']:.17g} | {sq['d50']:.17g} |
| d90 | {q['d90']:.17g} | {sq['d90']:.17g} |

Target mean 用 `sum(p_i*(low_i+high_i)/2)`，target 分位数用解析逆 CDF；sample 分位数采用 NumPy `method="linear"`，仅作诊断。**冻结 histogram 的连续均值约 2.069088 µm，与参考 2.51 µm 相差 {v['target_minus_reference_mean_um']:.9f} µm。** 2.51 的角色仅为 INFORMATIONAL_ONLY；没有为了消除此差异改概率或尾部。本 PASS 仅说明 sampler 重现了给定 histogram，不是对人工数字化或原始实验分布准确性的重新验证。

## 可复现性、文件及性能

两次独立 CLI 进程各抽 10000 个、seed=20260915，diameter 数组逐元素、float64 字节和完整 CSV 字节完全一致；临时 CSV 比对后删除，hash 回执保留在 [REPRODUCIBILITY_CHECK.json](validation/REPRODUCIBILITY_CHECK.json)。实测环境为 Python {platform.python_version()}、NumPy {v['environment']['numpy']} / {v['environment']['bit_generator']}、float64；不声称未测的 NumPy 版本或平台仍自动逐字节一致。

Histogram contract 不固定全局唯一 seed，`seed_policy=frozen_per_case`。每个 population 的相邻 `.metadata.json` 保存 seed、N、histogram SHA256、sampler version/source SHA256、环境和 CSV/数组 hash。CSV 的 source_bin_low/high 仅是逆 CDF 区间 provenance。V0 CLI 只允许 DEMO 或 VALIDATION role，拒绝隐式 seed 和覆盖已有 population。

单次 N=1000000、seed=20260917 的内存采样耗时 **{v['performance']['wall_seconds']:.9f} s**，约 {v['performance']['samples_per_second']:.3f} samples/s；包含新建 default_rng 与向量化逆 CDF，不含读取 histogram、写 CSV 或画图。未进行性能调优，百万样本仅用于计时及 hash 回执，没有保存或当作正式 case。

QA 图已生成并目视检查：

![Histogram QA](validation/SONOVUE_SAMPLER_QA.png)

![CDF QA](validation/SONOVUE_CDF_QA.png)

正式冻结合同：[SONOVUE_SAMPLER_CONTRACT_V0.json](contracts/SONOVUE_SAMPLER_CONTRACT_V0.json)。源码：[sonovue_sampler.py](src/sonovue_sampler.py)。完整数值核查：[SONOVUE_SAMPLER_VALIDATION.json](validation/SONOVUE_SAMPLER_VALIDATION.json)。完整性核查：在本目录运行 `sha256sum -c SHA256SUMS`（不自包含校验表本身）。

## 使用与阶段范围

本机已验证的解释器为 `/usr/bin/python3`。PATH 中的 Conda base Python 当前没有 NumPy；无需新安装，使用下面明确的解释器。

```bash
cd /home/lzy/projects/sonovue_size_distribution_v0
/usr/bin/python3 -B src/sonovue_sampler.py \\
  --histogram input/FROZEN_SONOVUE_HISTOGRAM.csv \\
  --n 1000 --seed 20260915 \\
  --output validation/ANOTHER_DEMO_1000.csv \\
  --role SAMPLER_DEMO_ONLY
```

API 可导入 `src.sonovue_sampler.sample_diameters(n, seed)`；默认读取随项目封存的 input 文件。`SonoVueDistribution` 同时提供 `cdf`、`pdf`、`inverse_cdf`。当前已有 demo 的 ROLE 为 SAMPLER_DEMO_ONLY，NOT_FORMAL_SIMULATION_POPULATION；正式 N/浓度将由未来 simulation case contract 决定。

MICROBUBBLE_TRANSPORT_MODEL、LAMMPS_PARTICLE_ENGINE、PALABOS_LAMMPS_COUPLING、MICROBUBBLE_WALL_MODEL、MICROBUBBLE_ADHESION 均 **PENDING**；RBC **OFF**。本轮未修改 HemoCell、Palabos、GPU Stage4、PBS/BSA numerics、RBC Stage1 或 vascular geometry。
'''
    (ROOT / 'SONOVUE_SAMPLER_REPORT.md').write_text(report)
    (ROOT / 'README.md').write_text('''# SonoVue size distribution V0

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
''')
    summary = {
        'SONOVUE_CONTINUOUS_SAMPLER_V0': 'PASS', 'PRIMARY_SOURCE': contract['primary_source'],
        'WEIGHTING': contract['weighting'], 'REPRESENTATION': contract['representation'],
        'CDF_TYPE': contract['cdf_type'], 'PDF_TYPE': contract['pdf_type'], 'SAMPLING_METHOD': contract['sampling_method'],
        'BIN_WIDTH_UM': v['bin_width_um'], 'SUPPORT_MIN_UM': v['support_min_um'], 'SUPPORT_MAX_UM': v['support_max_um'],
        'INPUT_HISTOGRAM_SHA256': v['histogram_sha256'], 'INPUT_PROBABILITY_SUM': v['INPUT_PROBABILITY_SUM'],
        'NORMALIZED_PROBABILITY_SUM': v['NORMALIZED_PROBABILITY_SUM'],
        'TARGET_MEAN_DIAMETER_UM': v['TARGET_MEAN_DIAMETER_UM'], 'VALIDATION_SAMPLE_MEAN_UM': v['VALIDATION_SAMPLE_MEAN_UM'],
        'TARGET_D10_UM': q['d10'], 'TARGET_D50_UM': q['d50'], 'TARGET_D90_UM': q['d90'],
        'VALIDATION_D10_UM': sq['d10'], 'VALIDATION_D50_UM': sq['d50'], 'VALIDATION_D90_UM': sq['d90'],
        'MAX_BIN_PROBABILITY_ERROR': v['MAX_BIN_PROBABILITY_ERROR'], 'RMSE_BIN_PROBABILITY': v['RMSE_BIN_PROBABILITY'],
        'MAX_CDF_DEVIATION': v['MAX_CDF_DEVIATION'], 'INVERSE_CDF_MAX_ERROR': v['INVERSE_CDF_MAX_ERROR'],
        'REPRODUCIBILITY_CHECK': 'PASS', 'SAMPLE_SUPPORT': 'PASS', 'PARAMETRIC_FIT': 'NO', 'SMOOTHING': 'NO',
        'REWEIGHT_TO_2P51': 'NO', 'DEMO_POPULATION_ROLE': 'SAMPLER_DEMO_ONLY',
        'SAMPLING_1E6_WALL_SECONDS': v['performance']['wall_seconds'],
        'MICROBUBBLE_TRANSPORT_MODEL': 'PENDING', 'LAMMPS_PARTICLE_ENGINE': 'PENDING', 'PALABOS_LAMMPS_COUPLING': 'PENDING',
        'MICROBUBBLE_WALL_MODEL': 'PENDING', 'MICROBUBBLE_ADHESION': 'PENDING', 'RBC': 'OFF',
        'NEXT_STAGE': 'MINIMAL_LAMMPS_PARTICLE_ENGINE', 'NEXT_STAGE_STARTED': 'NO', 'REPORT_DIR': str(ROOT),
    }
    write_json(ROOT / 'FINAL_SUMMARY.json', summary)
    (ROOT / 'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(f'{k} = {val}' for k,val in summary.items()) + '\n')
    for path in [contract_path, ROOT / 'contracts/VALIDATION_PLAN_V0.json', ROOT / 'input/FROZEN_SONOVUE_HISTOGRAM.csv']:
        path.chmod(0o444)
    files = sorted(p for p in ROOT.rglob('*') if p.is_file() and p.name != 'SHA256SUMS' and '__pycache__' not in p.parts)
    (ROOT / 'SHA256SUMS').write_text(''.join(sha(p) + '  ' + str(p.relative_to(ROOT)) + '\n' for p in files))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
