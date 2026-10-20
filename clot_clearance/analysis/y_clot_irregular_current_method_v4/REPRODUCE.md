# 原 Y 血块上的 V4 方法

本轮唯一正式轨迹 `20261020T041100Z`。从历史 Y 方法 v2 的完整原始参考几何迁移，V4 不规则点云、非均匀积分体积、变分力、不可逆软化、映射控制体积及表面载荷保留。原连续几何、原附着面、所有材料点 ID、相位 ID 与锚点坐标不变。边界材料点仅在初始扰动生成时保留原坐标；正式动力学中全部材料点可动，只有锚点固定。固定种子 20261019；保持 Y 案例原积分体积和质量，不使用长方体 V4 案例的质量。

运行位置 `/workspace/y_clot/project`；环境 `/workspace/y_clot/env.sh`。执行器和库仅用 1 rank / 1 线程，不依赖 GPU。依赖库为已冻结 V4 `verification/clot_irregular_lps_discretization_v4/20261019T221000Z/build/liblammps.so`。整个上游构建树不复制。

正式服务 `y_clot_irregular_v4_20261020T041100Z`；旋转视频服务 `y_clot_v4_rotation_20261020T041100Z`。均关闭自动重跑。`code/monitor.py` 只读查看当前进度和最后合法检查点。已有本轮输出时，不再次运行 `run_formal.sh`。

准备过程：`prepare_reference.py` → `prepare_boundary.py` → `check_elastic.py` → `checks.py` → `admit.py`。`SOURCE_FREEZE.json`、`V4_IMPLEMENTATION_FREEZE.json` 保存代码、材料、参考拓扑和原生库的哈希。正式过程 `run_case.py pilot`（到 2 ms）后，`run_case.py continue` 从同一检查点连续到 3.6 ms；它们组成同一条物理轨迹。数值不通过则停止并保留失败态，不自动调参或追加运行。

完成后 `terminal_audit.py` 只读恢复末态，不增加动力学步；`analyze_result.py` 计算存活原始图、碎块质量、形态、能量及与历史 Y 规则点阵的对照。历史案例读取既有轨迹，不重算。`render_results.py` 从实际接受态生成红色视频与载荷场图；`render_rotation.py` 固定物理时间 t=0，只让相机绕世界 Z 轴旋转一圈。

视频没有位移放大、节点插值、表面填洞或生成式图像。球半径由每点 Vi 换算，只是积分点标记。载荷图使用保存状态重放表面求积并核对真实外力，箭头只做统一比例显示。所有原始场保存在 `ACTUAL_LOAD_FIELD.npz`。

初版后处理后，`render_detail_stills.py` 收紧初态、大碎块阶段及载荷静态图的相机视野；原图与原 SHA 清单留在 `provenance/before_static_detail_refinement`，`STATIC_DETAIL_RENDERING.json` 核验全部载荷数组不变。之后再次执行 `finalize.py`、`validate_delivery.py`。动力学视频仍采用固定的全程视野，旋转视频不变。

科学范围：原几何为工程构造的三瓣壁附血块；并非病人重建。强弱材料尚未标定。1000 Pa 是指定当前表面切向滑移的机械载荷代理；本轮没有新解 CFD、声场、微泡 streaming、双向 FSI 或 5 s 输运。L1/L2 仍是部分当前表面近似。
