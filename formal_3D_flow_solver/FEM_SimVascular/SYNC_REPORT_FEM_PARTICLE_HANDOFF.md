# Stage FEM-FREEZE-SYNC 同步报告

## 从哪里同步？

唯一 source of truth：`/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular`。所有 accepted 资料来自 WSL；没有从 GPU 服务器重新下载，没有运行 CFD。代码/测试/配置逐文件复制并校验，原 WSL Git 历史和工作树保持原状。

## 同步到哪里？

GitHub：`OZAKI39/ulm-3d-vascular-model-generation`，分支 `sync/fem-simvascular-stage-q-particle-handoff-20260920`。WSL 独立克隆生成 commit，实际 remote 为 `https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`。

## main 有没有被修改？

**NO**。分支由 fetch 后最新 `origin/main` 的 `8a264014a6181d23d30f3ad4447931ce182e1ac4` 创建；main commit date `2026-08-23T18:51:21+02:00`。仅向新 sync 分支 push，无 force。差异限定为本 FEM 子目录及 root README 的短入口链接，原 main 代码不改写；最终审计见 `sync_metadata/main_history_audit.json`。

## 冻结 FEM 是哪一阶段？

Stage SV1.3Q，`ILU_REBUILD_POLICY_WINNER_FOUND`。RA 完整稳态 1571.003470 s，首次正式稳态 step70、安全停止 step71。这只是 frozen particle-development background-flow baseline，未扩大 scientific approval。

## 最终背景流是哪一个文件？

[frozen_reference/flow/steady_flow_stage_sv1_3q.vtu](frozen_reference/flow/steady_flow_stage_sv1_3q.vtu)，8997629 bytes，SHA256 `373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1`。

## Particle-0 需要哪些文件？

| 用途 | 相对此目录的路径 | SHA256 |
|---|---|---|
| 背景流 | [frozen_reference/flow/steady_flow_stage_sv1_3q.vtu](frozen_reference/flow/steady_flow_stage_sv1_3q.vtu) | `373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1` |
| 正式体网格 | [frozen_reference/SV_MESH/mesh-complete.mesh.vtu](frozen_reference/SV_MESH/mesh-complete.mesh.vtu) | `1a7a69dcba52f0475b2280775921e4dedcdfefa67b7965096d86cb9dba38bdf9` |
| WALL | [frozen_reference/SV_MESH/mesh-surfaces/WALL.vtp](frozen_reference/SV_MESH/mesh-surfaces/WALL.vtp) | `06c5d1d1975cf9470096d9b7e560e2173c9733c847d5bc364aab7355ef12a11b` |
| OUTLET_03 | [frozen_reference/SV_MESH/mesh-surfaces/OUTLET_03.vtp](frozen_reference/SV_MESH/mesh-surfaces/OUTLET_03.vtp) | `9e1d8d387eaff15ea3421c85051dd2f194d8adf57c8a839b4d894c5957e8546e` |
| OUTLET_01 | [frozen_reference/SV_MESH/mesh-surfaces/OUTLET_01.vtp](frozen_reference/SV_MESH/mesh-surfaces/OUTLET_01.vtp) | `1a051e6e2313c9be3985b16bd2a060fa111e8e8f8035682ab29eaa1961acdfed` |
| INLET | [frozen_reference/SV_MESH/mesh-surfaces/INLET.vtp](frozen_reference/SV_MESH/mesh-surfaces/INLET.vtp) | `fba5b444da15ebb7b881612acde9de1056395ff4c2903dd0b52e45e57b3fa9f5` |
| OUTLET_02 | [frozen_reference/SV_MESH/mesh-surfaces/OUTLET_02.vtp](frozen_reference/SV_MESH/mesh-surfaces/OUTLET_02.vtp) | `7dac9192aafbe935a595f4424c952e6fb6922e879b8e3821e3dc314f5b48cf33` |

同时读取 `frozen_reference/mesh_manifest.json`、`frozen_reference/boundary_manifest.json`、`frozen_reference/flow/flow_field_manifest.json` 和 `frozen_reference/baseline_summary.json`。单位、字段、法向、cell ordering 约定均在交接文档。Checkpoint 是额外可复现证据，Particle-0 不必读取它。

## 哪些内容没有同步？

PETSc/MPI/CUDA source/install/build trees、solver executable、CMake cache、编译中间文件、Python 环境、缓存、core、下载包、远端 home、SSH key/token/credential，以及巨大机器 inventory/重复中间运行输出。上游 svMultiPhysics 源码及其许可证保留在 vendor；950 个未下载历史示例/安装包的 LFS pointer 原文保存在 `upstream_lfs_pointer_inventory.json`，没有启用或下载 LFS。

完整候选扫描、复制映射及遗漏清单见 `sync_metadata/copied_from_wsl.json`、`sync_metadata/source_manifest.csv`、`sync_metadata/omitted_artifacts.json`。Stage Q full solver log、steady history、adaptive trace、stop request/ack、原始验收与环境摘要都保留。旧环境 inventory 用 hash pointer 替代，不将 65–82 MiB 的机器清单放入普通 Git。

## 有没有大文件被 pointer 替代？

必要 scientific artifact **没有**：最终 flow、正式 mesh、全部边界和最终 checkpoint 都是完整 Git 文件。历史机器清单、依赖 archive、重复中间输出的 WSL path/size/SHA 记录在 `omitted_artifacts.json`；取回时从其记录的 source 路径复制并重新核对 SHA，禁止把未取回文件声称为已验证。上游已有的 950 个未下载示例/安装包只有指针：其原文/OID/逻辑大小另列于 `sync_metadata/upstream_lfs_pointer_inventory.json`，这些从来不是正式冻结 scientific artifact。没有 zip 科学二进制绕过大小限制。文件大小检查见 `sync_metadata/file_size_audit.json`。

## 有哪些 validation 明确没有做？

Mesh convergence: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

Time-step sensitivity: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

CPU/GPU field-level equivalence: **DEFERRED BY USER DECISION**。

这些是 USER-ACCEPTED PROJECT DECISIONS，不是通过的验证。FEM performance tuning: frozen after Stage Q。Stage Q is the frozen particle-development background-flow baseline；不能解释为网格收敛、时间步收敛、CPU/GPU 科学等价或新的 FEM production approval。

本次执行 frozen artifact fresh-read、同 run config/time/hash 一致性、边界法向与流量复核、来源 manifest、无 secret、大文件及 main 历史审计；运行同步测试和原 Stage Q artifact tests，均不求解 CFD。测试记录见 `sync_metadata/test_summary.json`；原完整历史失败仍见 Stage Q 原 pytest 报告。

## 下一个聊天从哪里开始？

[PARTICLE_HANDOFF.md](PARTICLE_HANDOFF.md) → Particle-0。此次未写 FrozenFEMField、point locator、interpolation、RBC/MB dynamics、LAMMPS 或 resistance solver。

## 发布验收

本阶段最终状态由 `sync_metadata/stage_status.json` 记录；远端 push 与 GitHub 内容重新读取的证据在 `sync_metadata/remote_verification.json`。远端验证按已发布 commit 锚定，最终 commit 的 local/remote HEAD 一致性另由交付终端记录验证，避免把文件自己的 commit SHA 递归写入自身。
