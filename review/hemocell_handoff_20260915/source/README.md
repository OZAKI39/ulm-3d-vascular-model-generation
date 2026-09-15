# 当前实际源码映射

以下文件按已有 restore manifest、execution SHA256 和结果归档校验表选择。原件字节不改；各文件原路径、原哈希见 `../provenance/ORIGINAL_FILE_HASHES.tsv`，完整 bundle 映射见 `../SOURCE_MANIFEST.tsv`。

| 路径 | 已实际执行的职责 | 范围 / 依赖 |
|---|---|---|
| `gpu_stage4/source/vascularPoC.cpp` | 接受的 Stage4 短基准驱动 | 原 SHA256 `e23363ce5cedd9278e39bedec749ae90a38ff595a263d16e8e2e7a7dbb97e508` |
| `gpu_stage4/source/gpu_adapter.hpp` | Guo wall / velocity inlet / 三 pressure outlets、批量 halo、device safety/mass/flux、metadata residency | 相关逻辑已融合于真实 adapter；没有另造 pressure/Guo 文件 |
| `gpu_stage4/source/mpi_support.hpp` | MPI ownership / scalar reduction 与 flux 支持 | 与当前生产源一致 |
| `gpu_stage4/source/stage4_lbm_batch.hpp` | persistent descriptors + LBM batch dispatch | Stage4 accepted kernel 调度 |
| `gpu_stage4/source/stage2_profile.hpp` | 保留的低开销 profiling 范围 | 正常性能运行 profile=off |
| `gpu_stage4/scripts/` | 原实施、执行、评价与归档工具 | CMake 见同级 source；NVHPC26.5 原 cc120 |
| `rtx4090_restore/scripts/recipe/CMakeLists.txt` | 实际 cc89 恢复构建配方 | 同目录控制器/prepare_work 记录 Palabos 静态库和 TBB 路径 |
| `step3c_formal/source/` + `scripts/` | 正式长程控制、独立 CPU evaluator、在线监测、finalizer | source/formal_control.hpp 配合外部 evaluator_service；旧介质合同在 `../contracts/step3c/` |
| `new_medium/scripts/prepare_numerics.py` | 当前唯一 dt、出口 rho_LU 和 solver_parameters 生成源 | SHA256 `ca64f80c9be523a3466b20874a22842f0e1501b98a510d67f778f157984bfee3` |
| `new_medium/source/vascularPoC.cpp` | 当前 PBS/BSA smoke driver | SHA256 `8b71dd8d1d7fbc7e73b76de012de43f0c74ef05f54631caa7bc6947d2a2a2d65` |
| `new_medium/scripts/` + `inputs/` + `provenance/generate_driver.py` | 真实生成/构建/run wrapper/monitor/评价/打包实现 | 不把 provenance 中旧 generator 视为当前权威 |
| `rbc_stage1/source/` | 原生 reference RBC 网格工具及 HemoCell CPU library recipe | 仅 mesh construction，不是耦合 timestep binary |
| `rbc_stage1/scripts/` | mesh generation、确定性 pose search、独立 winding/distance audit、合同生成、finalize | 真实执行的静态工具；NumPy / VTK 已有环境 |

RBC runtime monitor/output policy 仅在 `../rbc_stage1/RBC_STAGE1_CONTRACT.json` 定义；**runtime 实现不存在**。`finalize_rbc_stage1.py` 核验的是本次静态构建/几何失败证据，不能用它声称 solver finalizer、MPI4、IBM 或 RBC safety 已验证。

## 补丁链与 upstream

GPU Palabos pinned commit 为 `4127697e90169bbef982295f1d1c933cf6e90caa`。原始 Stage4 全部变更在 `../patches/STAGE4_PATCH.diff`（从 Stage3-equivalent baseline 到 Stage4，含项目源和 4 个 native headers）。metadata 预置优化在 `STAGE2_MANAGED_METADATA.diff`，Stage4 adapter 已含该逻辑。**不要在已应用的生产导出树上重复打补丁。**

为了从原 GPU Palabos 导出恢复 native 头文件，本次还生成 `PALABOS_GPU_CUMULATIVE_STAGE4.diff`：仅含对比归档 pinned export 与已验证 Stage4 export 的 4 个不同 headers。完整共同文件 hash 对比和 patch 回放见 `../provenance/PATCH_REPLAY_VERIFICATION.json`。恢复 native 可以使用此累计补丁，然后放入已归档的最终项目源；不需要重新运行历史优化脚本。

旧 HemoCell + CPU Palabos 是另一条 upstream 链，不能与 GPU Palabos 混合：HemoCell `5a410848bd5c57d5ae1c171112e78eab4a82e650`，原 Palabos `05712164d940a42e06afdd705249912fa0c49f14` + `HEMOCELL_OFFICIAL_PALABOS.patch`。其原生 source/native 1512 个文件不复制到 GitHub；谱系证据在 `../rbc_stage1/provenance/SOURCE_LINEAGE_AUDIT.json`。

## 恢复注意

本次没有构建或运行任何代码。历史 CMake 和控制脚本保留实际 paths、只读来源核验和防重复运行逻辑，不是已完成移植的入口。要恢复，应在新的工作目录放入对应 upstream 和依赖、复制正确版本的 source/scripts，再根据原 source→work manifest 映射必要输入。new_medium 的 inputs 在本代码目录，生成合同在 `../new_medium/`；RBC 的合同/几何在 `../rbc_stage1/`。这种审查布局与原运行布局不同，不能通过删掉 provenance 断言来强行运行。

`multiplane_quadrature.tsv` 约39 MiB及原生体素数组未提交，可凭遗漏索引取回；control_volume_indices 与中小型合同已保存。原 STL/Step2 合同优先参考现有 `review_bundle/step3_review/`，必须核对身份而非重新几何生成。未来 RBC 运行仍需先通过当前失败的几何门并实现兼容壁面作用。
