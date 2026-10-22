# 本轮实际验证结果

以下 PASS 指对应的具体检查，不是整套流场物理精度通过。FAIL/BLOCKED 的旧科学结论完整保留。

| 检查 | 状态 | 实际范围 / 证据 |
|---|---|---|
| 冻结资料逐字节身份 | **PASS** | 两套场、全部旧群体、输入、历史证据与原生源码；[记录](../validation/FINAL_IMMUTABLE_DATA.json) |
| 原来源文件未改动 | **PASS** | 3439 个原路径的 SHA 复查；[记录](../docs/provenance/ORIGINAL_FILES_UNCHANGED.json) |
| 原 Git 工作树 | **PASS** | 相同参数核对 HEAD、状态数量和既有源码变更列表；[记录](../validation/ORIGINAL_WORKTREE_STATUS.json) |
| 几何、单位、两姿态 | **PASS** | 包含完整源表面坐标往返与速度大小保持；[记录](../validation/GEOMETRY_AND_POSES.json) |
| 原网格检查 | **PASS** | 节点、单元、边界、质量和封闭性；[记录](../validation/read_only_v1/mesh.json) |
| 冻结 CFD / WSS / checkpoint | **PASS** | 只读复核；没有新求解，不认证近壁精度；[记录](../validation/read_only_v1/cfd_wss.json) |
| 旧群体状态和账本 | **PASS** | 保留 141 观察截止、2 数值失败及 false 科学完整性状态；[记录](../validation/read_only_v1/population.json) |
| V2 已存误差分解复算 | **PASS** | 原科学结论仍为 MIXED_FIELD_ERROR；[记录](../validation/reference_reproduction/REVIEW.json) |
| 源码与数据独立性 | **PASS** | 隔离 Python，拒绝旧工作区读取，实际导入与数据载入；[记录](../validation/INDEPENDENCE.json) |
| 真实 BraVa 场 CPU/GPU 求值 | **PASS** | 2146 位置；无轨迹积分，不能证明局部物理准确；[记录](../validation/real_gpu_field/RESULT.json) |
| 网格软件加载 | **PASS** | SimVascular 内嵌 Python 导入 TetGen；没有生成新网格；[记录](../validation/MESHER_ENVIRONMENT.json) |
| 原生求解器构建和加载 | **PASS** | 实际编译、动态库和无输入启动；没有执行 CFD；[记录](../validation/NATIVE_ENVIRONMENT.json) |
| 当前报告和可视化链接 | **PASS** | 新报告、新索引、新 HTML；历史 HTML 的旧资源链接另列；[记录](../validation/CURRENT_LINKS.json) |
| 原 CPU/GPU 小型测试及工程保护测试 | **PASS** | {'total': 74, 'passed': 74, 'failed': 0, 'errors': 0, 'skipped': 0}；人工小场景，没有正式 BraVa 群体；[记录](../validation/small_tests/pytest.log) |
| 新正式 BraVa CFD / 新网格 | **NOT_RUN** | 本轮禁止启动；仅准备了计划；[记录](../validation/prepared_case_example/PLAN.json) |
| 新固定时间窗正式群体 | **NOT_RUN** | 协议与小测试已另列；没有新正式批次；[记录](../validation/prepared_population_example/PLAN.json) |
| V2 局部流场精度 | **FAIL** | 保留原科学判定 REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC；[记录](../diagnostics/reference_pipe/FIELD_VALIDATION.json) |
| P1 轨迹敏感性 / 收敛认证 | **FAIL** | 保留 P1_FIELD_TRAJECTORY_IMPACT_LARGE；并非本轮新增失败；[记录](../diagnostics/reference_pipe/P1_FIELD_PARTICLE_SENSITIVITY.json) |
| P2/P1 正式可用性 | **BLOCKED** | 原资源限制停止，且曲面压力检查仍有失败；本轮未恢复；[记录](../diagnostics/p2p1/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md) |
| V3 新网格和新场 | **NOT_RUN** | 保留提案身份；[记录](../diagnostics/reference_pipe/PROPOSED_V3_MESH_POLICY.json) |
| 可选旧 LAMMPS 接口 | **NOT_RUN** | 不是当前 BraVa 群体所用后端；未安装测试；[记录](../docs/NATIVE_ENVIRONMENT_ZH.md) |

完整测试使用本地实际 CUDA 设备；是否有跳过以 JUnit 数字为准。初次缺少依赖的收集失败、CPU 先行测试中的 GPU 跳过、安装超时及重试都保存在 `validation/logs/`，最终结果使用完整环境下的独立测试目录。

本轮没有完整 CFD 数值回归；原生二进制能启动，只是环境与接口检查。也没有把旧 1565 条结果当成新协议的大批量运行。

实际命令入口见 [README](../README.md)，新输出目录已经存在时会拒绝覆盖；复查请选择新目录。
