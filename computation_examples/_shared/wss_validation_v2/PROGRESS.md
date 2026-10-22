# WSS第二轮最终执行记录

## 最终范围和结果

1. 阶段一已完成：生产入口/材料读取修复，H0原始面片和显示回归差均0Pa，判别力测试及错误输入拒绝均实际运行。
2. 阶段二已完成：三档圆管实际CFD及最细dt/2。最细速度L2、中央压降、平均WSS误差分别0.292940%、−0.101695%、−2.834119%；达到原预登记目标。已完成验证不重算。
3. 阶段三已完成**两档网格敏感性分析**：原H0身份/稳态复核，中档CPU8/LU原标准从零求解至第51步、4079.429秒、118次线性校正无失败/遗漏/重试。独立生产WSS/日志/截面/控制体/几何检查通过。J1均值增加3.320%，J2 P95增加5.707%；最大内部截面偏差1.97302%Qin→1.11683%Qin；局部跳变及阈值小区域仍敏感，不能声明网格无关或确定离散误差。
4. 细档真实血管与阶段四压力敏感性均为**SKIPPED_BY_USER**。用户限制时间成本，细档未初始化/求解，O2两例未生成/排队/启动。细档网格/配置及取消记录保留，不当作已完成CFD。

## 取消操作和后台核查

细档`vessel_fine_LU16`及其他细档启动项已停止/禁用。压力计划的旧本地自动衔接进程已停止，活动顺序脚本已移除全部压力生成/注册/启动/等待逻辑，新增阶段取消保护。中档monitor PID617253从未因取消操作重启，正常完成退出。结果回传与独立分析已完成退出，没有等待取消阶段的流程。

证据：`stage3/vessel_fine/reports/user_cancellation.json`、`cancellation_verification.json`、`stage4/user_cancellation.json`、`remote_cancellation_verification.json`、`final_no_new_CFD_verification.json`。最终范围见`evidence/final_scope.json`；最初三档血管/细档初始化/压力计划保留为历史，已被用户最新指令替代。

## 失败和修复记录

- 首次圆管坐标保护过严：确认VTU精确float32序列化，仍用原double坐标求梯度；保留失败例，修改后已从零重跑并回归H0。
- 原网格GPU8及中档GPU1未恢复线性失败未作有效结果。已停止的CPU/GPU/MPI/ILU6/Schur性能实验均保留，不当作稳态验证。
- 原网格GPU1半步长一次旧预条件器失败后，原RHS恢复/零初值/新预条件器重试证据齐全；首次监视器FAIL记录保留，分类修正有8项实际/负测试。无容差放宽。
- CPU8/LU及CPU16/LU同第5步对照实际通过后采用CPU8/LU中档；细档CPU16计划未执行。MPI slots与实验分块配置首错保留。
- 停旧monitor未必终止独立MPI的问题已定位，按专用case cwd精确清理并加入信号处理；未影响其他项目或主中档。详情在`stage3/EXECUTION_UPDATES.md`。

## 交付证据

主报告`WSS_VALIDATION_REPORT_V2.md`、汇总CSV和五图PNG/PDF已生成；4份真实J1/J2场子集生产WSS恢复最大差0Pa。关键CSV包括两档区域、位置/阈值重合、邻片跳变、内部偏差及分流；压力响应没有伪造零数据。

有效运行、失败、停止、未执行、用户取消分别列入`CASE_STATUS.csv`；大文件本地/服务器路径在`ARTIFACT_INDEX.csv`。原保护文件270项最终核对见`evidence/protected_inputs_after.json`；精简ZIP及解压独立复核见`bundle/verification.json`。历史命令留在`evidence/COMMANDS_CHRONOLOGICAL_HISTORY.md`，当前有效复核入口为`COMMANDS.md`。
