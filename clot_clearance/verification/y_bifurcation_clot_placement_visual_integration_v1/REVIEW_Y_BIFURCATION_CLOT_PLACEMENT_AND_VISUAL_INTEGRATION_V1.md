# Y 分叉附壁血栓装配与可视化审查

结果：READY_FOR_HUMAN_REVIEW。本轮只构造一个 CURVED_MURAL_DOME_V1，完成实际 Y 网格导入、曲壁血栓、立方格材料点与三层壁锚点、新连接、零步初始化及背景 CFD 展示。没有启动含血栓 CFD、微泡或断裂动力学。CLOT_PRESENT_CFD_SOLVED = NOT_RUN。

1. 原始文件和身份

已通过 SSH 在服务器 `/workspace/mevo_acoustic_microbubbles` 逐级读取 config、geometry/mesh/steady/import STATE 及其实际引用文件，下载 15 项只读快照并逐项核对 SHA256。15 项也均与本地平行项目镜像相同。服务器所查项目根下没有当前 clot_clearance；新开发和零步检查在已有本地血栓工作树 `/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance` 执行，没有在服务器伪造血栓仓库或修改平行项目。

- 配置：`/workspace/mevo_acoustic_microbubbles/configs/stage4/bifurcation.yaml`；SHA256 `a1454e99da3ac1db8741d7b803e44d691d8cf0dc1566c3c4918171be0f201c1c`。
- 原四面体网格：`/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/mesh/20261002T141441884393Z/SV_MESH/mesh-complete.mesh.vtu`；SHA256 `74b7164dd981007557a86123b6854105e0302080c3e9b4685cf162aaacac10d3`。
- 原稳态场：`/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/steady/20261002T141810826990Z/run/8-procs/result_121.vtu`；SHA256 `131f7cb9b3e4ac9887bb0aa344d0eb2c0a2c8ac1d3b19fafb093ecd9a1045620`。
- WALL 与三个端口位于网格同目录 `mesh-surfaces/`；每个原文件的路径、大小、mtime、完整 hash 和只读快照路径见 SOURCE_ASSET_REGISTER.json / SOURCE_CONFIG.json。

ACTUAL_Y_ASSET_AUDIT.json 的 34 项手交信息均有 ACTUAL_FILE_VERIFIED 记录；配置来源和实际网格测量分别标注。母管长度 14 mm、半径 1.3 mm；两子管长度 12 mm、半径 1 mm、偏转 ±30°。真实体网格 26,712 节点、146,230 四面体。密度 1000 kg/m³、黏度 0.001 Pa·s、入口抛物线通量归一化和两个零参考自然牵引出口同时从实际 solver.xml 核对。原场第 121 步、dt=0.01 s，时间 1.21 s。此身份与 BraVa 四端口 18 mL/min 案例分开。

2. 放置位置、壁面拟合与间隙

中心线投影固定 (-6,0,0) mm，名义脚点 (-6,0,1.3) mm，实际 WALL 径向交点 (-6,0,1.297181237) mm。脚点修正 2.818763 µm，低于 dx=50 µm，不移动轴向位置。通过裁剪/线性细分实际 WALL 三角形获得连续足迹，每个足迹点沿局部径向向内使用指定 h=H(1−s²)^(3/2) 高度，未另造 Y 管。

足迹是一个连通面片；面积 0.991775 mm²。附着顶点对原壁面最大/中位距离为 0.000158817/0.000034924 µm。最大表面正向偏差 0.000158817 µm，小于预定 0.002 µm 数值容差；这是对同一离散表面的重合精度，不是生物壁面的纳米精度。材料中心离壁最小 28.450 µm，与体素/积分点离散有关，不能把中心距误报为连续实体悬空。

最终全部 1,364 材料点及 36,828 个体素中心/角点/面中点/边中点通过严格四面体归属；1,125 个锚点中心在原管腔外。表面、三端口与对侧壁的距离记录在 LUMEN_INTERSECTION_AND_CLEARANCE.json；入口最近顶点距离约 7.438 mm，未触及端口或分叉区域。锚点是固定壁代理点，未声称其名义体积构成真实血管壁。

3. 新形状与旧平面半球

| 指标 | 曲壁附壁隆起 | 旧平面半球 |
|---|---:|---:|
| 名义 a、b、H / 半球 R | 0.5625 mm | 0.5625 mm |
| 连续实体体积 | 0.192860 mm³（实际闭合网格） | 0.372757 mm³（2πR³/3） |
| LAMMPS 积分体积 | 0.170500 mm³ | 0.371375 mm³（历史保存结果） |
| 材料点数 | 1364 | 2971 |
| 附着形态 | 圆柱实际曲壁上的连续椭圆足迹 | 平面圆盘 |

新连续体积约为旧解析半球的 51.74%，不能称等体积变形。足迹轴向跨度 1.124855 mm、名义圆柱展开周向跨度 1.124604 mm；保存网格最大高度 0.562175 mm。外露面积 1.367675 mm²，总表面积 2.359450 mm²。材料中心 bbox、连续实体 bbox 与质心见 GEOMETRY_SUMMARY.json。

连续实体中心截面 x=−6 mm 面积 0.309141 mm²，除以同位置实际管腔截面 5.281726 mm² 得 5.8530%。这是几何截面占据率，不能当水力阻力或流量变化。另列材料柱的 dx² 占据率 5.2066%。

连续边缘平滑降至零；有限 dx 无法填满全部薄缘。保守掩码去除 104 个不满足壁面支持点条件的候选体素，另去除 1 个单邻居细须。积分体积相对连续网格相差 11.59%，没有进行分辨率收敛或擅自改 dx。医学实体视图显示生成材料点所用的同一连续参考网格，求解审查视图显示实际保留点和连接，二者差异明确保留。

4. LAMMPS 装配和零步读取

- 材料点 1,364；固定锚点 1,125，三层分别 358/390/377。
- BULK 57,346；INTERFACE 16,186；WALL 37,691；总连接 111,223。
- 三层按真实壁面外侧 0–50、50–100、100–150 µm 法向距离壳层选点，所有坐标保持全局 dx=50 µm 立方格；不移动格点贴墙。每点体积 dx³=1.25e−13 m³，每点质量 1.25e−10 kg。血栓总质量 1.705e-07 kg，排除锚点。
- δ=150.5 µm；无重复点、重复连接、孤立点；实际最大连接长度 150 µm。10 个候选 BULK 跨自由面连接被剔除；界面连接必须穿过实际 attachment patch 且不穿过 exposed surface。全部 WALL 连接另经精确线段/管壁相交检查，穿腔数为 0。
- 复用原 write_typed_data、LPS prelude、SI 内力密度到节点力转换与现有变体，未复制第二套求解器。变体二进制 SHA256 `7cba88869fb61531b4e5fe5bbcf3b129f80cc1d9c19c0a4761d55e69b67cf95b` 未改变。
- 最终 `read_data` + `run 0` 成功，2,489 个 atom、111,223 个 bond 与 NPZ 逐项对应，只有 step=0、time=0，坐标不变，速度/损伤/初始力为零。无时间积分 fix 或外加载荷。

开发记录保留：初步 VTK FindCell 曾接受一个略越壁的支持点，严格重心坐标检查拒绝后增强到每体素 27 点及 1 nm 保守裕量，最终去除该体素及一个细须点。先前草稿的 run 0 和图保存在 development_attempts，正式数据以本目录顶层 FINAL_STATUS 为准。总计两次零步初始化、零次动力学；没有多位置或多形态扫参。

5. 冻结 CFD 展示能说明什么

流线与箭头从原 VTU 速度数组生成，显示原管腔内的流向、速度分布，以及新血栓相对于原流场的位置。56 条显示流线中 4 条在第一次碰到新固体时截断，全部保存线段中点均在新固体之外。没有构造绕流、涡旋或假微泡。所有背景图与交互视图明确标记 UNOBSTRUCTED BACKGROUND CFD — NOT CLOTTED FLOW FIELD。

不能据此得到血栓造成的压降、分流变化、近血栓壁面剪切、表面准确 traction、微泡聚集或血栓破裂。background_flow.vtu 保留原无血栓场；新实体内部的原场只代表历史数值，不能作为现态流体。对可见流线做 mask 不会产生新的物理解。

6. 原流量和分流是否保留

从实际端口网格及原速度重新积分，入口约 1.58124363 mL/min；A/B 约 0.79056450/0.79067913 mL/min，仍为原近各半分流。全域原速度范围 [0.0, 0.010074017175615313] m/s，压力范围 [-0.0011369699976425888, 0.7060202132883139] Pa。原字段和导出背景字段的坐标、Velocity、Pressure 数组逐项相等。这些全部是原无血栓结果。

7. 展示和审查来源

四张 PNG/PDF、72 帧/6 s 镜头环绕 MP4/GIF、交互 HTML 和 ParaView VTM 均离线读取同一份保存坐标。仅在显示时统一 m→mm，各方向同尺度；没有单独放大或移动血栓。光滑效果仅插值表面显示法向，未移动实体点或 LAMMPS 点。局部视图从管腔侧显示隆起，并提供真实 y=0 剖面。

自动检查了 PNG 非空、PDF、MP4/GIF 帧、ParaView 数据重载、源哈希与图形数组对应；四张主图和动画四视角拼图人工/模型目视检查通过。HTML 已导出内嵌场景，未进行浏览器交互自动化；ParaView 数据已实际重新读入。首次交互导出缺插件，只在本任务 visual_dependencies 隔离安装 trame-pyvista / nest_asyncio2，未改变平行项目环境。

8. 未来载荷 provider 契约

SURFACE_NODES_AND_NORMALS.npz 给出 323 个初始外露面材料节点的全局坐标、外露法向、实际表面积权重 A_i、到真实壁面的距离、投影后的轴向切向、周向切向、壁到管腔方向、原局部流向及 geometry/mesh/field 身份。只将外露面三角形面积分配一次到最近外露材料节点，A_i 总和 1.36767506515e-06 m²；不加入附着面面积，不把全部体积节点当受载表面。薄缘节点可有两种拓扑邻近标签，面积集合仍互不重复。

provider 必须输出 load_source、surface_node_id、position_xyz_m、time_s/时间序列、source_geometry_id、source_flow_id、coordinate_frame、source_units、validity_domain，以及 force_xyz_N，或 traction_xyz_Pa 和 nodal_area_m2。后者逐节点用 F_i=t_i A_i 转换。未来沿真实曲壁切平面投影方向，不能继续用旧平底半球的原点径向。当前数据为初始映射，不支持断裂新增表面自动更新。

BACKGROUND_Y_FLOW = REQUIRES_CLOT_PRESENT_CFD_FOR_TRACTION；ANALYTIC_STREAMING 仅 schema；BUBBLE_STREAMING / CFD_STREAMING = NOT_YET。本轮没有把任何这些载荷施加给 LAMMPS。

9. 材料、接触和界面适用范围

保留 Liu 同源血栓未损伤小应变映射 K=33333.333333333336 Pa、G=940 Pa，工程密度 1000 kg/m³，bulk s00=0.15、interface s00=0.140624859375、α=0，以及已存在的 MODIFIED_1P0DX 接触变体。1.0dx 是有限工况下的工程接触选择，不是生物接触距离。界面是同文献猪皮剥离 proxy；曲壁法向、邻域和连接组织均变了，旧平面宏观强度不能自动迁移为本场景认证。

CURVED_INTERFACE_STRENGTH_VALIDATED = NOT_ESTABLISHED；材料映射、曲壁界面、接触及像素/网格展示都不构成临床、患者特异或实验精度验证。固定锚层不代表可变形血管壁或 FSI。

10. 真正的含血栓 traction 需要什么

下一阶段若获授权，应以原管腔减去血栓实体形成新的流体域，为暴露血栓表面建立独立 boundary tag，重建相容体网格及近表面解析度；移除附着区域的流体空间并核对端口和流固边界拓扑。可把原入口流量与出口牵引作为清楚声明的初始 BC，但必须重新求解新的含血栓速度/压力，再由总应力 σ=−pI+μ(∇u+∇uᵀ) 及统一法向约定映射表面牵引，并检查表面积、总力、单位和局部误差。若血栓会变形或断裂，另需定义几何更新与耦合策略。

原 CFD 的求解域包括现在的固体区域，且未在新 exposed surface 上施加流固无滑移边界，因此不能把旧场压力/梯度直接命名为新血栓上的真实 traction。以上仅列出必要后续，不在本轮执行。

交付与停止

- results 目录 OPEN_RESULTS.html：四图、视频、交互 HTML 和 ParaView 入口。
- verification 顶层 FINAL_STATUS.json、输入 data/deck、run0 日志、全部几何/流场/source 检查及 evidence.zip。
- 历史 1775 个文件、14064505446 bytes 哈希保持；15 项服务器和本地镜像原资产保持。没有 commit/push 或未完成本轮作业。
- 现阶段主用途：EXPLORATORY_GEOMETRY_AND_BACKGROUND_FLOW_INTEGRATION。停止并等待人工审查，不自动执行 CFD、bubble 或 fracture。

实现参考：[LAMMPS BPM/PERI](https://docs.lammps.org/bond_bpm_peri.html)、[read_data](https://docs.lammps.org/read_data.html)、[PyVista UnstructuredGrid](https://docs.pyvista.org/api/core/_autosummary/pyvista.unstructuredgrid)。官方 BPM/PERI 的等间距立方格、节点体积和读取规则用于输入装配；工程接触变体的身份单独保存。
