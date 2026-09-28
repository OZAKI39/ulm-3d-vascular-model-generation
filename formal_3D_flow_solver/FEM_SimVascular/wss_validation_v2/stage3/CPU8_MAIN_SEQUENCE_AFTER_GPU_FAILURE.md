# 中档GPU1失败后恢复已登记CPU8序列

中档GPU1实际启动后，第一轮 fresh 线性求解在300次迭代 DIVERGED_BREAKDOWN；第二次逻辑求解的 fresh retry 亦在100次迭代失败，明确 FRESH_RETRY_FAILED。此错误不能按成功恢复解释，未获得可接受CFD场。停止本次专用进程组，原输入、原日志、监控及失败状态保存在 stage3/failed_attempts/vessel_medium_gpu1/，没有覆盖原始H0。旧GPU中/细等待程序已停用。

根据此前已完成的CPU8/GPU1同场校验，恢复中/细网格生成后、GPU选择前已经登记并存档的CPU8配置。恢复后每个 input_hashes.json 逐项核对，完全等于 reports/unexecuted_cpu8_inputs 中的原登记版本；网格、XML、dt、PDE、黏度、边界和收敛容差均未改变。两档均从零初值实际求解。

vessel_medium_cpu8_retry 等待原网格CPU8同条件完整复算通过，然后启动中档；vessel_fine_cpu8 等待中档通过。当前基线仍为已合格原H0，CPU8独立复算核对执行一致性；中/细不再采用GPU1。此决定有实际失败证据，不是根据WSS高低选择结果。GPU后端故障尚未唯一定位；它不证明网格或物理模型错误。

为避免未恢复线性失败继续浪费计算，验证监视器新增早停：发现 FRESH_RETRY_FAILED 或已经接受到NS日志中的未收敛警告，即停止本次专用求解进程组并记录FAIL；正常成功重试不触发。原始日志不删除。所有成功场仍需完整稳态、代数残差、输入身份及后处理检查。
