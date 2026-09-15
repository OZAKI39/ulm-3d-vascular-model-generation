# 当前状态（2026-09-15）

本表汇总截至本次交接的状态；科学失败不会因 GitHub 归档成功而升级。

| STEP | STATUS | EVIDENCE | NEXT ACTION |
|---|---|---|---|
| STEP1 geometry | PASS | [冻结几何](../../review_bundle/step3_review/step1/geometry_contract/README.md) | 保持单位、STL 与几何身份 |
| STEP2 voxelization / ports | PASS | [Step2 报告](../../review_bundle/step3_review/step2/STEP2_REPORT.md) | 保持体素和四端口映射 |
| PURE_FLUID_BC | PASS | [物理 BC](contracts/step3c/physical_bc_contract.json)；[旧 Step3 报告](../../review_bundle/step3_review/step3/STEP3_REPORT.md) | 保持 Guo / 入口 / 三出口数学 |
| HUMAN_PARAVIEW_REVIEW | PASS | [用户当前声明与历史状态范围](provenance/STATUS_RECONCILIATION.md) | 不扩展为 RBC 人工检查通过 |
| STEP3C / PURE_FLUID_BASELINE | PASS | [旧介质正式长验证](pure_fluid/STEP3C_FORMAL_GPU_REPORT.md) | 保留原 Qtarget 与倍率；新介质未长期校准 |
| GPU_STAGE4 | PASS | [Stage4 报告](gpu_stage4/original_rtx5090/GPU_STAGE4_REPORT.md) | 保留生产 candidate；停止自动纯流体优化 |
| NEW_VAST4090_PLATFORM | PASS | [恢复报告](gpu_stage4/rtx4090_restore/VAST4090_STAGE4_RESTORE_REPORT.md) | 使用实际 cc89 构建来源；不复用未知 binary |
| PBS_BSA_NUMERICS | PASS | [smoke 报告](new_medium/PURE_FLUID_NEW_MEDIUM_SMOKE_REPORT.md)、[合同](new_medium/NEW_MEDIUM_NUMERICS_CONTRACT.json) | 实验测量/最终介质合同 PENDING |
| RBC_NATIVE_LIBRARY_AND_MESH_BUILD | PASS | [构建 provenance](rbc_stage1/RBC_STAGE1_BUILD_PROVENANCE.json) | 构建成功不代表耦合 runtime 验证 |
| RBC_STAGE1 | FAIL_GEOMETRY_GATE | [Stage1 报告](rbc_stage1/RBC_STAGE1_REPORT.md)、[几何报告](rbc_stage1/RBC_GEOMETRY_FIT_REPORT.md) | 解决入口加载与初态几何；RBC timesteps=0 |
| RBC_WALL_INTERACTION | ABSENT | [源码审计](rbc_stage1/provenance/WALL_INTERACTION_AUDIT.json) | 建立与 Guo STL 兼容的壁面作用 |
| HUMAN_RBC_REVIEW | PENDING | [几何诊断](rbc_stage1/geometry/README.md) | 人工查看同坐标系的 lumen 与失败姿态 |
| RBC_STAGE2 | PENDING | [当前机器状态](CURRENT_STATE.json) | 完成兼容性阶段并重新冻结/验证 Stage1 后再议 |
| RBC_GPU | PENDING | [Stage1 最终状态](rbc_stage1/FINAL_SUMMARY_LOCAL.json) | 尚无 GPU IBM / RBC mechanics 证据 |
| HCT_5_10_20_PERCENT | PENDING | [开发定义与研究次序](NEW_CHAT_CONTEXT.md) | 未启动；不可从静态网格推断群体正确性 |
| MICROBUBBLE | PENDING | [最终研究目标](NEW_CHAT_CONTEXT.md) | 当前不加入 microbubble |

**NEXT_STAGE = RBC_GEOMETRY_AND_WALL_COMPATIBILITY_STAGE**：upstream loading / inlet injection、Guo-compatible RBC wall interaction，以及 nominal 50 vs actual 45.046 µm³ 定义审阅。
