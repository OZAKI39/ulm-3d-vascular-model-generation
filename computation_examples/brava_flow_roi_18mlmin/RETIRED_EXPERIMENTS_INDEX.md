# Retired experiments index

只在明确历史审计时读取。下列最终结论在本次清理提交归并；表中commit是清理前实现恢复点，不表示该commit包含所有未跟踪原始数据。

| 路线 | 最终状态 | 最终报告 | 本次清理前commit | 一行原因 |
|---|---|---|---|---|
| P1 brute-force refinement | RETIRED_NO_PRODUCTION_WINNER | [FINAL_VALIDATION_REPORT_ZH.md](mesh_particle_validation_v2/FINAL_VALIDATION_REPORT_ZH.md) | `5622c1a9549c` | 局部零速缺陷可消除，空间输运仍未收敛 |
| Local P2/P1 | LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED | [FINAL_P2_NEARWALL_REVIEW_ZH.md](p2p1_nearwall_cross_validation/FINAL_P2_NEARWALL_REVIEW_ZH.md) | `5622c1a9549c` | A区近壁稳定性门槛未通过 |
| BC-B / frozen Jacobian plan | INCONCLUSIVE / PLAN_NOT_RUN | [FINAL_BOUNDARY_DIAGNOSIS_ZH.md](boundary_transfer_diagnosis/FINAL_BOUNDARY_DIAGNOSIS_ZH.md) | `5622c1a9549c` | NS线性求解受限，边界gate不可评估；Jacobian计划退休 |
| TetGen scalar metric / NW1 | RETIRED_RESOLUTION_GATE | [NEARWALL_TETRA_CALIBRATION_ZH.md](production_workflow/NEARWALL_TETRA_CALIBRATION_ZH.md) | `5622c1a9549c` | 合法NW1高度50.16µm未达20–30µm；baseline保留 |
| MMG anisotropic / runtime | CLOSED_FOR_NOW | [MMG_RUNTIME_DIAGNOSIS_ZH.md](production_workflow/MMG_RUNTIME_DIAGNOSIS_ZH.md) | `5622c1a9549c` | A1、D0/D1/D2在有界预算内无最终网格 |
| Omega_h | RETIRED_RUNTIME_LIMIT | [TRAJECTORY_MESH_RESULT_ZH.md](production_workflow/TRAJECTORY_MESH_RESULT_ZH.md) | `5622c1a9549c` | 目标metric步之前原生修复膨胀，1800s停止 |
| NASA/refine | PRECHECK_FAILED_DEFERRED | [TRAJECTORY_MESH_RESULT_ZH.md](production_workflow/TRAJECTORY_MESH_RESULT_ZH.md) | `5622c1a9549c` | 当前CAD关联退化；R1未运行 |
| MeshSim | UNAVAILABLE_KERNEL_INTERFACE | [TRAJECTORY_MESH_RESULT_ZH.md](production_workflow/TRAJECTORY_MESH_RESULT_ZH.md) | `5622c1a9549c` | 现有kernel不可实例化；license状态未确立 |

恢复本次删除的tracked实现：`git show 5622c1a9549c19306b28eb5479c87600d76997c0:PATH`。必要日志片段和tool/config SHA见`RETIRED_EXPERIMENT_EVIDENCE.json`；不恢复后台任务。

更早已退休的实现用`git log --all -- PATH`定位；首轮清理前源码基线为`11dfa45763df55f51dffdd702c1e23ed65751c2e`。本次不声称5622c1a包含更早已删除的实现。

## 2026-10-04 同步清理

本轮清理V1/cap旧绘图、旧CPU monitor、中间GPU step10–70、旧CPU瞬态场/重复检查点、per-second监督trace及重复PDF/SVG；最终V2和GPU第71步保留。保留原BREP/INLET/BOUNDARY_TOPOLOGY/RAW_MESH与CPU step10 checkpoint，因为未改动的当前生成/QA/prepare仍依赖它们；这些不构成继续旧路线的授权。当前三图保留PNG。

GPU完整日志2.82MB与195KB/5秒resource samples均属当前identity直接依赖，予以保留，另加GPU_STEADY_RUN_EVIDENCE.json。旧CPU简报见 production_workflow/manufactured_pipe/R1_GMSH_V2/cfd/reports/CPU_TIMEOUT_SUMMARY.json。历史报告对已删除产物的路径均只作 retired 证据；不得把它们当当前运行入口。

本轮tracked实现恢复点：本地独立算例Git `91ebc02162e7b0c774cf649b5669bcc5c8026e54`；GitHub同步仓库 `0b808fb6c90e27674d79078c49358ced2e9258ff`。二者不是同一Git历史。未跟踪原始数据不能声称可由Git恢复；原服务器只读保留，清理清单记录路径与SHA。第三方完整checkout从GitHub活动树删除，当前GPU wrapper/RA header/policy/build身份与host-sync patch保留于gpu_solver_fix。
