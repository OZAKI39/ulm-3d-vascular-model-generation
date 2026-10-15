# 当前clot原生软化开发入口

当前：lps_native_softening_small_nonaffine_crack_validation_v1 / 20261014T214309Z（用户明确取消4 GiB RSS停止条件后的同工况新运行）。

- 状态 **STOPPED_WALL_CLOCK_BUDGET**；实际到2.4ms，9600接受/0损伤拒绝；计划2.5ms完成=False。停止原因["AssertionError('WALL_OR_OUTPUT_BUDGET')"]。
- 3080点、144448原始参考键、3524预切、700×500×1100µm试样及所有原始m/r0/V/taper不变。上下各140点仅规定z、2800内部全自由；每侧120µm、2ms半余弦+0.5ms保持不变。
- 仅取消RSS停止门槛，无替代内存上限；原1800s墙钟/2GiB输出及所有数值门槛保留。复用SHA730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90隔离库，无C++新编译、无P1新物理运行。
- 旧50µs仅FP64而完整restart在0，故按用户授权从原始预切t=0重新推进。前201标量行/11节点帧与旧0–50µs逐值一致，无前缀拼接。旧资源失败及P1漏检历史保留。
- 本次自然软化16548，新增完全失效14342；首软化0.0007475000000000181s、首完全失效0.000799750000000025s；区域分类NOTCH_LOCALIZED_DAMAGE。存活跨面3312，分量节点数[3080]。
- 最大接受Δd=0.078523184524，最大归一化余额0.15671857463%；全部3080原始ID保留。数值准入、损伤分类与新生裂面认证分开。
- CRACK_AREA_UNRESOLVED，G_SEPARATION_PROTOCOL/GC_EFFECTIVE均UNAVAILABLE；宏观强度、全局dt收敛及V3生产仍未认证。
- 末态仅50µm无桥接开口带，未达到100µm持续前缘推进准则；严格horizon无桥接子集为空。裂尖先软化，但末态大量失效位于端部/其他自由表面，不能把分类名理解成连续裂纹传播已认证。
- 最新完整原生检查点时间2.4ms：/workspace/y_clot/project/verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T214309Z/full/checkpoints/accepted_slot_0.restart及其.state.npz，已有SHA。若下一轮补足剩余保持段，应恢复原生历史与Python能量账本，不能直接把现有fresh-t0入口称为续算；本轮未追加恢复运行。
- 主段墙钟1826.010s，CPU1749.435s，峰值RSS11.24126GiB；1 rank×1线程，GPU未用。资源详情见RUN_RESOURCES.json。
- 先读verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T214309Z/FINAL_STATUS.json、34项REVIEW和NUMERICAL_VALIDITY.json；代码analysis/small_crack_no_rss_v1/。已有输出禁止覆盖或自动追加工况。
- 六组PNG/SVG：results/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T214309Z/OPEN_RESULTS.html；全部真实点阵/事件/曲线，无医学外壳或位移放大。
- WSL实际/home/lzy/projects/clot_rupture.md已更新并读回，before/after/diff/SHA完整保存；①–⑨历史与旧G0标称值保持。
- 341点700µs历史成功基线保留：lps_native_softening_domain_retention_and_failure_safe_output_v1 / 20261014T150825Z。生产库、env.sh及其他历史数据不变。
- 下一步依据本次具体事件与现有数值停止原因审查，不自动放宽门槛、不改材料/几何/载荷、不运行V3/CFD/微泡/AI。前入口备份在本run的ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md。
