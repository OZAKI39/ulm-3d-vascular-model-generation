# Particle-0 → Particle-8 研究路线

FEM 已冻结。以下仅是研究路线，本次没有实现任何一阶段；统一采用第一版单向耦合 FEM → particles。

- **Particle-0 — Frozen FEM sampling**：定义冻结场读取、四面体点定位、插值和速度梯度接口，验证单位、边界和域外行为。
- **Particle-1 — single MB**：研究单个微泡在冻结背景流中的运动与基本力学，建立可检查的单粒子案例。
- **Particle-2 — single rigid RBC / orientation**：建立单个刚性 RBC 的平动、转动与取向模型，不延用旧 reduced-order RBC 模型作为正式动力学。
- **Particle-3 — wall interaction**：研究最近壁面、法向与粒子壁面相互作用，明确接触与穿透处理。
- **Particle-4 — particle contact**：研究粒子之间的接触与约束，并验证多体接触一致性。
- **Particle-5 — near-field hydrodynamic resistance**：研究近场流体阻力及其与接触模型的配合，避免重复计入作用。
- **Particle-6 — LAMMPS bridge**：研究与 LAMMPS 的数据及积分接口，明确职责、单位和同步方式；本阶段不安装它。
- **Particle-7 — inlet/outlet + continuous injection**：实现阶段将研究入口通量、持续注入、出口离域及可复现随机性。
- **Particle-8 — full RBC + MB suspension**：组合经验证的模块研究 RBC 与 MB 悬浮体系，报告适用范围和数值限制。
