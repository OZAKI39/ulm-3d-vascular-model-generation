# 半步长成功重试的监视器复核

GPU1 原网格半步长实际 CFD 已于第51步退出，exit_code=0，耗时1812.276秒。首次监视器只要发现任意 DIVERGED 标记便判 FAIL；它没有区分被丢弃的失败线性尝试与随后接受的校正。

原求解器 preconditioner_reuse.h 28–72行在复用预条件器失败时复制回原 RHS、VecEqual 检查、核对零初值，然后重建预条件器并重试。实际 logical=3 attempt=0 在500次迭代 DIVERGED_BREAKDOWN；日志随后明确 original_rhs_restored=1 initial_guess_zero=1，attempt=1 在554次迭代收敛，healthy=1 recovery=STALE_ILU_RECOVERED，二者之间没有接受 NS 校正。独立解析120次最终接受的线性校正均成功，未解析行0，非线性末相对残差6.542e-15，最大true residual/停止阈值0.9992582。原失败尝试计数仍为1，不能写成全部尝试零失败。

本轮只修改独立验证监视器 run_solver.py，新增 log_acceptance.py；没有修改求解器、PETSc 容差、CFD配置或解。七项真实/变异日志测试核对：正常圆管和成功重试接受；真实GPU8未恢复失败、缺少RHS恢复证明、缺少零初值证明、重试未收敛、提前接受校正均拒绝。见 data/monitor_recovery_classification_tests.json。

重新判定还要求原输入哈希、输出manifest哈希、exit_code、连续稳态间隔、末态字段和独立日志质量全部通过。首次判定保存在 reports/execution_original_monitor.json；新的 execution.json 内保留 divergence_markers 和验证后的重试详情；reports/monitor_reassessment.json 给出原日志/分类器/独立质量检查哈希。服务器和本地分别运行同一复核脚本，CFD没有重跑。独立GPU8失败例继续拒绝，未并入有效结果。
