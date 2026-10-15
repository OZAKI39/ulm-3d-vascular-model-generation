# 当前clot原生软化开发入口

本轮：lps_native_softening_event_resolved_integration_v1 / 20261013T235945Z（2026-10-14 Europe/Zurich）。

- **局部损伤事件PASS，完整700µs FAILED_ATOM_LOSS。最后接受687.75µs，下一688µs试算固定域丢点。不得继续或称完整协议通过。**
- 原bond1813的0.11363469跳跃用真实VV两段0.125µs分辨，最初两段Δd=0.05736736/0.05626720；局部15接受/4拒绝，正式2755接受/4损伤拒绝，另1次未完成原生候选。
- 拒绝后x/v/F/κ/d/U/D/计数/time、实际存活键和完整special排除均精确恢复。正式35个检查点，未每步写restart；预注册最小0.03125µs，实际最小0.125µs。
- 最大接受Δd=0.0978676344，接受区间最大相对能量余额0.372426823%。全局时间收敛仍未认证。
- 卸载到650µs完成，末保持仅37.75/50µs。末态4918完全失效、7129存活；原图连通分量192/147/1/1。Gc不可用。
- 丢点诊断：原盒子fff、各轴±250µm。自由节点73最后y=249.953987866µm、vy=0.400969997m/s、F=0；下一漂移预测y=250.054230365µm。原生日志实报340vs341。预测不是接受态。
- 因原生异常后gather/restart失败，最终restart、全键稀疏事件流、精确在线误差极值和正式峰值RAM未落盘。逐步标量及18个完整FP64关键态仍在；禁止用NPZ制造restart。
- 先读verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z/FINAL_STATUS.json、24项REVIEW、NATIVE_ATOM_LOSS_DIAGNOSIS.json、EVENT_INTEGRATOR_DESIGN.md。
- 活跃实现analysis/event_resolved_v1/；冻结执行代码SHA在EVENT_INTEGRATOR_FREEZE.json。run_protocol.py保留真实失败路径，禁止在封存目录重跑。
- 图表results/lps_native_softening_event_resolved_integration_v1/20261013T235945Z/OPEN_RESULTS.html（五组PNG/SVG）。
- 研究库verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z/build/liblammps.so，便携副本binary/liblammps.so；SHA 08ccf24110e2cf2d07dbc8c0f31c0f368368b19b5427e6e5de9cc108a5ae115d。
- 仅fix_move.cpp限定六equal规定运动变dt，新增event_state_audit.cpp只读排除表。bond/协调fix/pair/nve/Verlet源码不变；输入G/K/ρ/dx/δ/s0/sc/α、m/V/r0/taper和半余弦加载均不变。
- 动态接触最大重叠仍0.246649196dx；未作接触刚度标定，物理精度未认证。
- 生产环境source /workspace/y_clot/env.sh保持；生产库SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a。旧归档102文件在本轮入口更新前全部SHA一致，旧manifest未改；before入口保存在本轮目录。
- WSL实际文档已更新并读回，凭据WSL_DOCUMENT_UPDATE.json及wsl_document/；所有旧文本和①–⑨保留。
- 下一轮优先LPS_NATIVE_SOFTENING_DOMAIN_RETENTION_AND_FAILURE_SAFE_OUTPUT_V1，先审查自由分离点的计算域保留与异常证据持久化；本轮不自动执行。
- P4非均匀小裂纹NOT_ADMITTED/NOT_RUN；未运行V3、CFD、微泡或AI，也未部署生产或更改材料。
