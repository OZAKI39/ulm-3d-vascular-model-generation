# Vast RTX4090 Stage4 恢复核查

结果：**FAIL_STOPPED**。新平台纯流体基线：**FAIL**。

归档重新检查 15344 个文件，工作副本 5654 个文件；终检源码/冻结输入差异 1。只适配工作路径、工具链/链接路径和实测 RTX4090 cc89。原 Stage4 C++、已应用 patch 的 Palabos 导出源码、数学/数值合同及运行参数字节保持不变。原 RTX5090 binary 仅记录来源，未执行。

NVHPC26.5 / CUDA13.2 从冻结 NVIDIA 官方地址下载，SHA256 与原下载 provenance 相同；安装证据在 toolchain_env 和 provenance。未修改主机驱动或全局 CUDA。OpenMPI 仅安装缺失依赖；正式 GPU smoke 前完成无求解器的 MPI 单核绑定检查。GPU activity 由官方 cavity 的 CUDA kernel 原始 trace 证明；Stage4 性能运行不 trace、不加同步、不优化。

| 验证 | 结果 | 说明 |
|---|---|---|
| 官方 cavity3d | NOT_RUN | 最多200步；kernel UNVERIFIED |
| Stage4 200 | NOT_RUN | 从零；RTX参照0/1/10/100/200、CPU参照0/100/200 |
| Stage4 1000 | NOT_RUN | 从零；0/100/200/1000全场及每步24组flux/质量/安全 |
| Stage4 5000 | NOT_RUN | 从零；0/100/200/1000/5000全场，执行原增长判据 |

每个 horizon 都对冻结的 CPU MPI12 和 RTX5090 Stage4 两组历史真实结果逐行比较。完整检查压缩CSV、各量最大误差、检查点数据和未修改比较器已归档，阈值未改变。NaN/Inf累计：NOT_RUN；趋势：UNVERIFIED。归一化速度、密度、总质量/CV质量、全部多平面流量与独立 MPI ownership、安全证据均由原比较器审查。

| 性能范围 | 初始化 s | 求解 steps/s | 端到端 s |
|---|---:|---:|---:|
| 1000（计时100–1000） | None | None | None |
| 5000（计时500–5000） | None | None | None |

持续求解速率相对 RTX5090：None；相对 CPU MPI12：None。分母为用户指定历史参考 125.880863 和17.692465 steps/s。每个 horizon 仅一次新测量，不宣称重复实验统计可靠性。初始化及端到端包含所需输出；性能不作为数值正确性门槛。

5000步内存判据：NOT_RUN；性能稳定性：NOT_RUN。具体 steady thirds、256MiB/512MiB阈值、首末1000步窗口及原Stage4 5% horizon下降判据见 RTX4090_PERFORMANCE.json。GPU利用率为一秒采样、按已完成步数筛选，不能用来宣称kernel占比。未运行10000步，所以不声称通过原5000/10000跨horizon重复稳定性测试。

执行状态：START。错误（如有）：AssertionError()。

远端 work：`/workspace/hemocell_restore/work/stage4_restore_20260915_150938`；compact：`/workspace/hemocell_restore/results/stage4_restore_20260915_150938`。原始全场二进制保留在远端各run，compact带全文件SHA256清单和完整逐项比较证据；较大的CSV仅在compact副本无损gzip，原run文件未改。官方原始CUDA trace随compact保留。WSL下载后单独验证manifest，远端summary此时尚未宣称本地传输PASS。

旧 incomplete upload：`/workspace/hemocell_restore/archive/vast_snapshot`，3564896256 bytes（allocated），未删除。正式历史基线和只读archive未改；没有RBC、PBS/BSA、microbubble、Stage5或long run。下一步仅作建议：REVIEW_FAILURE_BEFORE_ANY_NEW_RUN，本任务不会执行。
