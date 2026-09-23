# Particle-6.5 近场正则化文献审计

核对日期：2026-09-21。以下按原论文/出版记录核对，使用转述。文献先例、跨领域证据与用户批准的项目参数分别记录，不把它们混成一次实验测量。

| 来源 | 可核对的文献信息 | 本轮采用的事实 | 角色 |
|---|---|---|---|
| A | Christopher Ness, *Simulating dense, rate-independent suspension rheology using LAMMPS*, Computational Particle Mechanics **10**, 2031–2037 (2023), 2023-06-05 在线发表；[DOI](https://doi.org/10.1007/s40571-023-00605-x) | 第 3 节式 (3) 后，润滑的典型低间距限制为 O(10⁻³ a)，外截断为 0.05 a。 | DIRECT_SUSPENSION_SIMULATION_PRECEDENT |
| B | Keliu Wu 等, *Wettability effect on nanoconfined water flow*, PNAS **114**(13), 3358–3363 (2017)；[DOI](https://doi.org/10.1073/pnas.1612608114)、[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5380095/) | 引言讨论直径大于约 1.6 nm 的纳米孔中，连续介质与分子描述结合的可能适用性；边界条件、界面作用与表观黏度不能忽略。 | CROSS_DOMAIN_CONTINUUM_SCALE_EVIDENCE |
| C | Wentian Zheng 等, *Experimental observation of liquid–solid transition of nanoconfined water at ambient temperature*, Nature Materials **25**, 495–501 (2026)；2026-01-12 在线发表；[DOI](https://doi.org/10.1038/s41563-025-02456-8)、[PubMed 41526516](https://pubmed.ncbi.nlm.nih.gov/41526516/) | 摘要报告 hBN/亲水金刚石之间的受限水在约 1.6 nm 开始明显结构和扩散变化，低于约 1 nm 观察到完整结晶。 | CROSS_DOMAIN_NANOCONFINEMENT_EVIDENCE |
| D | D. Marsh, *Water adsorption isotherms and hydration forces for lysolipids and diacyl phospholipids*, Biophysical Journal **55**(6), 1093–1100 (1989)；[DOI](https://doi.org/10.1016/S0006-3495(89)82906-6)、[原文，PMID 2765647](https://pmc.ncbi.nlm.nih.gov/articles/PMC1330575/) | 摘要中二酰基磷脂酰胆碱的水合力衰减长度约为 0.3 nm（流动层状相）或 0.2 nm（凝胶相）。 | LIPID_INTERFACE_SCALE_CONTEXT |
| E | Juan Tu、Jingfeng Guan、Yuanyuan Qiu、Thomas J. Matula, *Estimating the shell parameters of SonoVue microbubbles using light scattering*, JASA **126**(6), 2954–2962 (2009)；[DOI](https://doi.org/10.1121/1.3242346)、[原文，PMCID PMC2803720](https://pmc.ncbi.nlm.nih.gov/articles/PMC2803720/) | 单气泡壳参数估计部分采用恒定 4 nm 壳厚假设。作者没有单独分析该厚度参数；不能写成本文直接测得的普适壳厚。 | SONOVUE_INTERFACE_SCALE_CONTEXT |

PNAS 论文的 [2019 年更正](https://pmc.ncbi.nlm.nih.gov/articles/PMC6717266/) 补充了先前会议版本和图 1–4 的来源。本轮已核对这条更正，未将其解释为新的接触长度测量。来源 C 的 DOI 含 2025，但正式在线和卷期记录属于 2026。

## 直接模型先例的适用边界（DIRECT MODEL PRECEDENT）

来源 A 是致密颗粒悬浮液模拟先例。其惯性、摩擦和接触力实现不是本项目采用的动力学；本项目仍复用原过阻尼阻力平衡与无摩擦速度约束，LAMMPS 仍只管理状态、邻居、输出和重启。数值尺度先例不等于 SonoVue–小鼠血管体系的直接标定。

## 跨领域证据的适用边界（CROSS-DOMAIN PHYSICAL EVIDENCE）

来源 B 的尺度是纳米孔直径，来源 C 使用特定固体表面，来源 D 是脂质体系的水合力衰减长度。几何、液体环境和界面材料均不能与本项目血管薄膜一一对应；这些结果只说明纳米尺度的体相连续介质解释需要谨慎。我们不引入这些论文的表面力、滑移、黏度函数或结晶模型。

## 本项目明确批准的选择（PROJECT MODEL CHOICE）

`NEAR_FIELD_REGULARIZATION_V1` 冻结 χ_full = 0.01、χ_off = 0.05、2 nm 分子尺度下限和 10⁻³ a_ref 悬浮液尺度，并取两者最大值作为 h_lower。2 nm 是 continuum handoff model choice，**不是直接测得的 SonoVue–小鼠血管 cutoff 或真实固体接触距离**。χ_full 和 C1 cubic smoothstep 属于本项目截断耦合选择；两球的 a_ref = min(a_i,a_j) 是保守参考长度选择，阻力系数仍用原 R_eff。

4 nm 仅作壳层背景，**不增加到 sampled MB radius**。本项目采样半径对应既有外径合同。糖萼模型和 RBC 润滑仍未冻结；2 nm 不代表糖萼厚度或几何偏置。1.5 nm 与 3 nm 只做敏感性比较，科学可接受性由用户审核。
