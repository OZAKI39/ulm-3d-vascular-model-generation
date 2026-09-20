# SonoVue 尺寸分布：GitHub 交接

本目录归档已经通过验证的 SonoVue 连续经验 CDF、直径采样器及尺寸分布图。
本分支从仓库 main 的 `8a264014a6181d23d30f3ad4447931ce182e1ac4` 建立；新增内容仅限本目录。

## 图表与入口

- [尺寸直方图：冻结输入与 100000 样本对照](validation/SONOVUE_SAMPLER_QA.png)
- [连续 CDF 与样本经验 CDF](validation/SONOVUE_CDF_QA.png)
- [采样源码](src/sonovue_sampler.py)
- [冻结 histogram](input/FROZEN_SONOVUE_HISTOGRAM.csv)
- [冻结 sampler contract](contracts/SONOVUE_SAMPLER_CONTRACT_V0.json)
- [中文科学报告](SONOVUE_SAMPLER_REPORT.md)
- [原阶段使用说明](README.md)

## 内容与科学范围

25 个原有文件（7448788 bytes）原样同步，包含代码、8 项单元测试、
冻结输入/合同、两张 QA 图、CDF 表、1000 个 demo 样本、100000 个验证样本、
相应 metadata、验证结果和原 SHA256SUMS。

分布为 fresh SonoVue、count-normalized / number-weighted；
原始 0.1 µm bins 的概率保持不变，采用 bin 内均匀密度及分段线性 CDF，
通过逆 CDF 生成连续直径。support 为 0.75–5.25 µm，来自冻结 CSV。
没有平滑、参数拟合、合并 bin、改变尾部或为匹配 2.51 µm 重新加权。
目标连续均值约 2.069088 µm；采样器 PASS 不代表重新验证人工数字化准确性。

这些 population 是 DEMO / VALIDATION，不是正式实验浓度定义。
本目录报告中的 PENDING / NEXT_STAGE 是原尺寸分布阶段的历史状态，
不代表整个微泡项目的最新开发状态。
本次同步没有运行流体、粒子输运、壁面或 RBC 模拟。

## 同步与可复现性

- 原 25 个文件与本地已验证目录逐字节一致。
- 原 SHA256SUMS 中的 24 个条目全部通过；校验表自身由同步 manifest 记录。
- 在新的仓库路径运行原有 8 项单元测试，全部通过，确认默认输入的相对路径有效。
- 新增本交接说明及 [GITHUB_SYNC_MANIFEST.json](provenance/GITHUB_SYNC_MANIFEST.json)。
- 原报告中 WSL 绝对路径作为历史 provenance 保留；采样核心默认输入相对模块目录定位。
- 原 validation/finalize 脚本保留科研记录，可能包含历史环境路径且会写输出；
  不应直接在冻结归档上重跑它们。

在已安装 NumPy 的 Python 环境，从本目录执行：

```bash
python3 -B -m unittest discover -s tests -v
sha256sum -c SHA256SUMS
```

API 与 CLI 用法见原 README 和科学报告。复现时记录每个 case 的 seed、
N、histogram SHA256 和 sampler version；不宣称未测试平台/NumPy版本的逐字节一致性。
