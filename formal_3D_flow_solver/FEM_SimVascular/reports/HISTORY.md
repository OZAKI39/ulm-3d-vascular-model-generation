# FEM 历史与证据范围

当前背景流唯一入口是 Stage SV1.3Q 的 frozen_reference，不能用旧流场替代。

- [SV1](sv1/REPORT.md) / [SV1.1](sv1_1/REPORT.md)：真实 SimVascular meshing 与线性求解稳定化。旧短算例质量守恒失败仍是失败证据。
- [SV1.2](sv1_2/REPORT.md)：完整 transient-to-steady 验证。
- [SV1.3](sv1_3/REPORT.md)：CPU early-stop；保留历史 production 身份，不在本阶段运行 CPU。
- [G](sv1_3g/REPORT.md)、[H](sv1_3h/REPORT.md)、[J](sv1_3j/REPORT.md)：GPU runtime、CUDA/PETSc 旧版兼容与矩阵测试；部分路线失败，被后续升级取代。
- [L](sv1_3l/REPORT.md)：修复缺失的 MPI Fortran predefined datatypes，继续暴露 CUDA ghost vector 兼容问题。
- [M](sv1_3m/REPORT.md)：旧 PETSc CUDA ghost backport，继续暴露 teardown/lifecycle 问题。
- [N](sv1_3n/REPORT.md)：现代 PETSc 3.25.5 + CUDA13.2，完成 lifecycle repair 和真实 GPU flow bring-up；历史 G/H/J/L/M 失败不因该结果被改写成成功。
- [O](sv1_3o/REPORT.md)：GPU 简单调优及 native final-output-on-stop。
- [P](sv1_3p/REPORT.md)：ASM2/ILU2、同时间步复用，保留 fallback。
- [Q](sv1_3q/REPORT.md)：跨时间步 reuse 比较及最终 RA 完整稳态，作为冻结的粒子开发背景流。

关键原报告、现有测试与小型证据保留；巨型机器 inventory、build 和重复中间输出不迁入 Git。遗漏路径与 SHA 在 `../sync_metadata/omitted_artifacts.json`。原 Stage Q 完整 suite 的 8 项历史失败见 `../PORTABILITY_NOTES.md`；本次未重跑 CFD 或将旧失败标成 PASS。
