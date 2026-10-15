# 当前clot原生软化开发入口

本轮：lps_native_softening_domain_retention_and_failure_safe_output_v1 / 20261014T150825Z（2026-10-14 Europe/Zurich）。

- **DOMAIN_RETENTION_PASS、FULL_700US_PASS、ENERGY_PASS、TOPOLOGY_PASS全部通过。** 唯一正式路径从原始无损伤t=0到700µs；未拼接旧前缀。
- 固定fff数值域由各轴±250µm扩大到±400µm，不remap、不加墙、不改周期性。依据旧真实范围、最大速度和余下时间规划。盒子不是物理血管壁。
- 全程341个原始ID保留；最小域余量145.134130µm。节点73在688µs实际y=250.054230µm，700µs为254.865870326µm，确认旧域丢点原因。
- 2804接受/4损伤拒绝，8个0.125µs接受步、2796个0.25µs步；冻结最小0.03125µs、实际损伤门槛0.10及原回滚算法保持。最大Δd=0.09786763438279739；最大相对能量余额0.372426823%。
- 正常原生末态gather核验341点、7110存活键与独立d<1原图端点集合完全一致；4937完全失效，4个实际分量192/147/1/1。节点73和105为孤立点。尚非已认证真实血栓碎块。
- W/U/K/D/Uc=4.575844829448/0.525683244103/0.674093509273/3.272865781378/0.126299310063nJ，Jc=0，R=-23.097015368980pJ。沿用完整变分能量账本和累计归一化。
- 旧0–687.75µs共2756标量时刻/18完整态比较相容；每步dt、失效数、接触数一致。旧逐步x/v缺失，只比较实际可用数据；未重跑旧域。
- 初始旧/新域物理状态逐值一致；683µs旧检查点上的一次异常输出测试通过，候选/接受态分开，失败处理不调用native gather或write_restart。
- 每步刷新FP64节点/事件/接触对流，周期两槽原生restart，沿用损伤检查点。最终有效原生restart与FP64配对见full/LATEST_NATIVE_CHECKPOINT.json。
- 先读verification/lps_native_softening_domain_retention_and_failure_safe_output_v1/20261014T150825Z/FINAL_STATUS.json、22项REVIEW及DOMAIN_PROTOCOL_FREEZE.json。新活跃驱动analysis/domain_retention_v1/，封存目录禁止重跑覆盖。
- 五组PNG/SVG：results/lps_native_softening_domain_retention_and_failure_safe_output_v1/20261014T150825Z/OPEN_RESULTS.html。真实点阵显示，不生成连续外壳。
- 复用旧隔离库verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z/build/liblammps.so，SHA 08ccf24110e2cf2d07dbc8c0f31c0f368368b19b5427e6e5de9cc108a5ae115d；本轮无C++改动、无重新编译。生产库和env.sh不变。
- 主计算墙钟54.473s、CPU35.475s、峰值659.688MiB；1 MPI×1线程，无GPU。完整资源与输出规模见RUN_RESOURCES.json/PAYLOAD_MANIFEST.json。
- 1.0dx接触保持，最大重叠0.246649196dx≈12.33246µm，物理精度尚未认证。全局时间收敛、宏观强度、材料Gc和完整V3生产就绪仍未认证。
- WSL实际/home/lzy/projects/clot_rupture.md已更新并读回；before/after/diff/SHA见WSL_DOCUMENT_UPDATE.json。全部历史内容及G0标称值保留。
- 下一轮唯一建议LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1；小型非均匀裂纹可研究，本轮NOT_RUN。未运行V3、CFD、微泡或AI，交付后停止等审查。
- 上一入口完整备份：verification/lps_native_softening_domain_retention_and_failure_safe_output_v1/20261014T150825Z/ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md；所有历史科研结果保留。
