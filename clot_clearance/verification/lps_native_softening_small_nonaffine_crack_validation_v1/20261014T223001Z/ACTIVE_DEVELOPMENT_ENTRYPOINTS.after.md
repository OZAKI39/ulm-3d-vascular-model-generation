# 当前clot原生软化开发入口

当前：lps_native_softening_small_nonaffine_crack_validation_v1 / 20261014T223001Z（取消1800秒墙钟预算后，从原生2.4ms检查点续算）。

- **FULL_PROTOCOL_NUMERICAL_PASS**：原定2.5ms时程已完整完成，累计10000接受步、0损伤拒绝；本次仅追加400步，无第二个物理工况、无t=0重跑。
- 3080点、144448原始参考键、3524预切、700×500×1100µm试样及完整m/r0/V/taper不变。上下各140点仅规定z；每侧120µm、2ms半余弦+0.5ms保持不变。
- 墙钟和RSS上限均为null，没有替代上限；2GiB输出门槛和全部损伤、能量、运动学准则保留。复用SHA730e11651aee6533b4fb7d38c8341e06ad853bf3397c4908fc6ef905d4e3ac90隔离库，无C++新编译或P1新物理运行。
- 恢复完整原生x/v/κ/d、原始ID/m及端点图，累计W/D/Jc、接触功、尺度、四项数值缺陷和区域D连续，能量R跳变量=0。旧0–2.4ms标量、节点流、损伤和接触流保持精确前缀。
- 本次末态自然软化16616、自然完全失效14392；首软化747.5µs、首完全失效799.75µs，均位于裂尖侧表面，原始bond ID68923/77580。裂尖软化领先加载端91.25µs。
- 最终存活跨面3312（长键3266），分量点数[3076, 1, 1, 1, 1]；主连通体仍有3076点，4个小分量是单节点，不能称为已形成两个宏观分离体。
- 最大接受Δd=0.078523184524、最大归一化余额0.157192777384%；原生末态端点图与独立原始图核对通过，全部3080原始ID保留。
- NOTCH_LOCALIZED_DAMAGE说明真实裂尖初始局部损伤；末态开口候选带仅50µm，未达到100µm持续前缘推进准则。后期加载端和其他自由表面大量耗能，不能由分类名认证连续裂纹传播。
- CRACK_AREA_UNRESOLVED，G_SEPARATION_PROTOCOL/GC_EFFECTIVE均UNAVAILABLE。网格代理0.025mm²、严格horizon无桥接子集0；不以总D除以代理面积。末态K/U=1.22225774，明显动态；最大接触重叠0.154408631dx，未认证接触物理精度。
- 最新完整原生检查点2.5ms：/workspace/y_clot/project/verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T223001Z/full/checkpoints/accepted_slot_0.restart及其.state.npz。原定窗口已结束，不自动延长物理时程。
- 本次续算墙钟68.625s、CPU65.197s、峰值RSS319.527MiB；原路径累计墙钟1894.635s、峰值RSS11.241257GiB。1rank×1线程，无GPU。
- 先读verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T223001Z/FINAL_STATUS.json、34项REVIEW及CONTINUATION_CHAIN_AUDIT.json；代码analysis/small_crack_resume_no_wall_v1/，新运行先检查冻结身份，禁止覆盖本归档。
- 六组PNG/SVG：results/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T223001Z/OPEN_RESULTS.html；使用实际点阵/事件/曲线，无医学外壳或位移放大。
- WSL实际/home/lzy/projects/clot_rupture.md已更新并读回，before/after/diff/SHA保存；原①–⑨、旧G0标称值及全部历史文字保留。
- 先前50µs内存停止、2.4ms墙钟停止历史保留；341点700µs成功基线20261014T150825Z保留。生产库、env.sh及历史物理数据不变；仅当前入口在保存旧字节后推进。
- 当前不认证宏观强度、全局dt收敛、材料Gc或V3生产；不自动改材料/几何/载荷、运行V3/CFD/微泡/AI。已完成本轮交付，等待下一轮科学审查。
