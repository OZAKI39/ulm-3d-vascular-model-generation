# 1. Executive summary

本轮为 **DISCOVERY ONLY**，事实源仅为当前 WSL 的 `/home/lzy/projects/ulm_3D_vascular`。HemoCell 仅观察结构和已有 Git 状态。未迁移、复制几何、运行源工程脚本、测试、solver、构建、安装、提交或 push。

**CONFIRMED_FROM_CODE / CONFIRMED_FROM_CONFIG / CONFIRMED_FROM_DATA**：当前主线是 s1–s5 分段流程，保存几何、局部/全局拓扑、显式端口标签和旧 APES 数值证据。**INFERRED**：适合先提取只读拓扑检查，再评审单位/坐标/边界契约；不能直接迁旧 CFD 配置。

事实标签：`CONFIRMED_FROM_CODE`=实际代码；`CONFIRMED_FROM_CONFIG`=配置声明而非实验真值；`CONFIRMED_FROM_DATA`=当前文件内容/头/元数据；`INFERRED`=推断/建议；`UNVERIFIED`=证据不足。单位未知使用 `UNIT_UNVERIFIED`，坐标未知使用 `COORDINATE_SYSTEM_UNVERIFIED`。本报告的迁移类别和风险等级均为 INFERRED。

范围：根目录深度 2 起步；发现重要目录后再定向读取。文件系统统计覆盖全树；数据扫描仅格式/大小/关联信息，共 18,296 条数据路径，未加载 28 GiB 数据体。静态解析 155 个主项目 Python 模块；完整读取 7 份配置；只检查代表性 TIFF/NPZ/H5/NIfTI 头和一个 1.64 MB 带标签 VTP；67 个采集 XML 与 ZIP 中央目录只读。PDF 全文、巨型压缩包内容、PKL 反序列化及第三方源码全面审阅未执行。

# 2. Source project identity

- CONFIRMED_FROM_DATA：path/realpath 均为 `/home/lzy/projects/ulm_3D_vascular`，与 `/home/lzy/projects/hemocell_starter` 独立。
- remote：`https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`
- branch：`codex/cfd-wall-force-numerics-validated-sync-20260830`
- HEAD：`fd21a850a16d0ba17ef3d521123badd53864a3a2`
- 普通 `git status --short` 因缺少 `git-lfs` 失败，原始错误完整保存于 `source_identity_raw.json`。没有安装 LFS 或修改 Git 配置。
- 只读替代：NUL 格式 status 临时禁用 LFS filter 后有 1620 项，335 项经 index LFS 指针的大小+SHA256 核对为已水合相同内容；剩余 1285 项：1114 修改、12 删除、159 未跟踪。完整列表在 `source_git_worktree_assessment.json`，不能把旧工程称为 clean。
- 这些是既有修改，本轮不 reset/restore/clean。旧入口删除与 s1–s5 新文件并存；结论以文件当前内容为准，不以 GitHub/HEAD 内容替代。
- Ultraliser HEAD `3e4b0eee685adbf513e40720a68fd92e66a34b44`，429 行 tracked status；LBPM HEAD `6d686d354e5b8140841d3601e4c8c0e4e4b77e48`，545 行 tracked status。仅保存身份，不判断这些差异均为算法修改。见 `nested_repository_identity.json`。
- 适用 `/home/lzy/.codex/AGENTS.md` 为空；源工程未发现额外 AGENTS.md。

# 3. Limited project map

CONFIRMED_FROM_DATA：总占用约 28 GiB，37,896 文件、3,046 目录；其中 vessel_model 文件内容约 22.71 GB，outputs 约 3.86 GB，.git 约 2.46 GB（这些是文件逻辑大小，不等同 du 占用）。

| 目录/模块 | 实际角色 |
|---|---|
| s1–s5 / configs / utils | 主入口、严格配置、几何/ROI/旧 CFD 模块 |
| vessel_model | 9 个数据集/归档/运行时混合资产族，不能整目录搬运 |
| outputs | 19 个全局 SWC 批次、19 个采样批次、3 个模型目录、12 个表面目录、42 个 CFD 目录 |
| test_data | ROI003274 最小 fixture，而非完整原始数据 |
| Ultraliser | 几何外部工具源码及 build-wsl |
| external_reference/LBPM | 另一个 LBM 参考仓库；不在当前 s1–s5 正式调用链 |
| docs / references / outputs/documentation | 说明、论文、教学图与历史证据，部分陈旧 |
| .git / .codex_tmp / caches / tmp | 不迁移的工程元数据与缓存 |


# 4. Main entrypoints

| 入口 | 实际调用 | 功能/限制 | 证据类别 |
|---|---|---|---|
| s1_swc_roi_generate.py:148 | load_swc_roi_yaml_config → run_rodent_vasculature_pipeline → run_sampling_from_rodent_run | 原始 SWC/TIFF/Mask → analysis SWC/图 → ROI NPZ、CSV | CONFIRMED_FROM_CODE |
| s2_swc_stl_model_generate.py:70 | load_swc_stl_yaml_config → resolve_sampling_run/select ROI → run_ultraliser_reconstruction | 保存 ROI → SWC/H5 → Ultraliser → STL/VTP/QC | CONFIRMED_FROM_CODE |
| s3_cfd_1D_data_preprocess.py:42 | load_cfd_preprocess_config → run_cfd_preprocess → solve_global_flow / transfer_all_boundaries | 全局图+ROI → 一维压力/流量及假定端口 | CONFIRMED_FROM_CODE |
| s4_cfd_surface_prepare.py:39 | load_surface_prepare_config → run_vmtk_surface_prepare → local_cut / VMTK / tagged caps | 局部切面+延长+重网格+标签及边界侧车 | CONFIRMED_FROM_CODE |
| s5_cfd_flow_solve.py:17 | load_cfd_flow_config → run_cfd_flow → promotion replay | FRESH_STEADY 分支实际抛授权错误；replay 可解析旧输出并可选 smoke | CONFIRMED_FROM_CODE |
| main.py:170 | run_pipeline / run_hierarchical_graph_pipeline | 旧 STL 清理→体素→中心线图；默认 unprocess_stl 路径不存在 | CONFIRMED_FROM_CODE |
| doc_visualize_v1/v2/v3.py | main / generate_visualizations / render_* | 保存结果可视化；v3 明确依赖旧稳态 VTU/物理截面字段 | CONFIRMED_FROM_CODE |


代码逐函数行号、imports、调用目标见 `code_structure_index.json`。入口脚本只被读取，连 `--help` 也未执行。C++ 重建入口位于 `Ultraliser/apps/ultraVessMorpho2Mesh/VessMorpho2Mesh.cpp:133`；它是外部几何程序，不是当前血流求解器。

# 5. Configuration chain

CONFIRMED_FROM_CODE：s1 → `utils/swc_roi_yaml_config.py:load_swc_roi_yaml_config`；s2 → `utils/cfd_lumen/model_yaml_config.py:113`；s3/s4/s5 分别调用各子包 `config.py` 的严格 YAML loader。大部分相对路径按 PROJECT_ROOT 解释，s2 的用户 YAML 参数先按 cwd 解析。源码使用 `yaml.safe_load`，校验字段/类型/互斥 ROI selector。

CONFIRMED_FROM_CONFIG：s2 和 s3 指向 `outputs/sampling/20260825_133201_radius_plus_structure_k5`，anchor=3274；s3 的 model_run/rodent_run 为 null，可按兼容性自动解析；s4 和 s5 固定特定旧输出路径。模型目录最新并不等于正式 CFD 输入。`configuration_inventory.json` 保存所有键值、SHA256 和具体路径是否存在。

INFERRED：顶层 `cfd_lumen_config.yaml` 属于旧通用配置入口参考，s2 实际读取完整 `configs/swc_stl_model_generate.yaml`；不能把两个配置当一份。VMTK 配置中的 D:/anaconda3 路径与 ../external/vmtk 在当前 WSL 不成立。当前配置 FRESH_STEADY 在 `utils/cfd_flow/pipeline.py:346` 的实际分支抛出 `CFD_FLOW_FRESH_STEADY_REQUIRES_EXPLICIT_COMPUTE_AUTHORIZATION`；不能声称默认 CLI 可直接跑完整新流场。

# 6. Input data

CONFIRMED_FROM_DATA：`data_file_inventory.json` 对每个匹配路径记录格式、字节、推测角色、代码/配置关联、单位/坐标证据和迁移价值；`data_family_profiles.json` 给出族级证据，族内未读文件不会自动继承代表样本单位。

- fMOST 活动样本：`.../raw_data/analysis_data/analysis_data/{images,mask,swc}/fMOST_0_5_6_0_0_6_0001_02_01.*`。image/mask 为 192³ TIFF，SWC 有七列；没有长度单位头。原图与 total_vascular_data 同名图像 SHA 相同，但 mask 和 SWC 不同，不可按名字去重。
- `test_data`：analysis NPZ 有 7419 节点；ROI NPZ 有 109 节点、108 边、3 CUT_PORT、1 TRUE_TERMINAL，保留 local/global IDs 和坐标/半径数组。它是测试输入，不是独立实测流场。
- 当前成功模型 H5：points=(113,4)，structure=(5,2)，connectivity=(4,2)，属性明确 coordinate_unit=um、第四列 diameter_um、radius_scale=0.91。
- 实际格式还包括 STL/VTP/VTU/VTI、CSV/JSON/YAML、TIFF/PNG/JPG、NIfTI、NPZ、H5、MAT、SWC、PKL、GraphML、旧 Lua/LSB/RES；泛用 NPY 示例出现在打包运行时，不要误认活动血管输入。OBJ 来自 Ultraliser 示例。
- `.msh/.ply/.dcm` 和未压缩传统 `.vtk` 在当前格式扫描中 NOT_FOUND；ZIP 内成员另记，不能当已解压输入。AneuX/Brain_arteries/IntrA 等存在压缩归档，中央目录见 `source_dataset_archives.json`，本轮未解压。
- 代表 NIfTI 头 `xyzt_units=0`、qform=0，虽有 sform/pixdim，仍不能猜 mm 或 um。Schmid PKL 有网络和 RBC_trajectories 文件，未反序列化，单位及实验/模拟来源 UNVERIFIED。

# 7. Geometry generation pipeline

CONFIRMED_FROM_CODE：原始 analysis SWC + TIFF/Mask（辅助） → xyz×spacing、保留原始半径 → 单连通 analysis SWC/层级图 → ROI 裁切及 CUT_PORT/TRUE_TERMINAL → ROIRecord/NPZ → canonical SWC + diameter H5（feed radius×0.91） → Ultraliser → lumen_surface_um.STL/VTP + lumen_surface_m.STL → 一维端口传递 → VMTK 局部切面/延长/重网格/封帽 → CellEntityIds/port_id VTP + 四端口 STL/CSV/JSON → 旧 APES 分区与笛卡尔 mesh。

原始 mask 是辅助核查/可视化数据，不会自动重写 canonical SWC 半径和父子边；s2 依赖保存的 ROIRecord，不是直接从任意 TIFF 重建。`utils/cfd_lumen/ultraliser_backend.py:248` 将 H5 radius feed 乘 0.91 并以直径写第四列；`ultraliser_qc.py:40` 导出 um 与 m 副本。ROI 构建保留全局物理位置，不把 local node ID 当局部坐标。

实际存在：lumen surface、STL/VTP surface mesh、centerline/SWC/GraphML、mask/TIFF、voxel/教学 VTI、ROI、branch labels、显式入口/出口 patch。旧 volume mesh 是 Seeder 笛卡尔 LSB 和导出的六面体 VTU；不是 HemoCell 可直接加载的通用四面体网格。原始未封帽、修补与封帽版本均有历史记录，不能仅取“最新 STL”。

逐项存在性、完整路径与用途限制见 `geometry_feature_inventory.json`。计算 ROI 已确认；独立实验 target region 为 UNVERIFIED。CONFIRMED_FROM_CODE / CONFIRMED_FROM_DATA：`docs/CFD_SURFACE_PREPARE.md` 描述旧环形延长/平滑算法，与当前 s4 的 VMTK 调用链不同，归为历史参考。选中历史表面是 cap-only recovery 产物，复用了已有 open surface；不能声称已由本轮当前代码完整重跑。

# 8. Output pipeline

CONFIRMED_FROM_CODE / CONFIRMED_FROM_DATA：各阶段按输出批次写 source config、manifest、geometry/data、QC、figures/report。s1 写 rodent_vasculature 和 sampling；s2 写 input/SWC/H5 与 geometry/STL/VTP、qc、report；s3 写 global_1d、roi/boundary_conditions.json 和 port planes；s4 写 geometry、boundaries、bc、vmtk、qc；s5 的 replay 解析已有 Base/LSB，输出稳态物理场 VTU、metrics CSV、QC JSON 和 production_review.html。

历史 `production_tau1_base_promotion_anchor003274_20260902_013637/qc/run_summary.json` 的状态是 `CFD_FLOW_PRODUCTION_TAU1_INTEGRATION_AND_VISUAL_REGRESSION_PASS`，execution_mode=`VALIDATED_BASE_PROMOTION_REPLAY`。这是读取到的历史记录，**不是本轮验证通过或新鲜求解**。docs/CFD_PRODUCTION_INTEGRATION_CONTEXT.md 的“尚未 promotion”描述已经落后于当前代码/数据。

CONFIRMED_FROM_DATA：当前 `docs/CFD_FLOW.md` 的 scientific scope 记录 Base/Coarse 稳态已接受，但仅为双网格分辨率敏感性；Fine 稳态和正式三网格 GCI 未完成，WSS 为 DEFERRED。本轮未复验，也不据此声称网格无关或壁面剪应力准确。该文档仍使用已删除的旧 CLI 名称，运行入口以当前 s5 源码为准；证据见 `historical_validation_scope.json`。

# 9. Units

| quantity | value | unit | meaning | confidence / source |
|---|---|---|---|---|
| 原始 fMOST SWC 坐标 | 原始三列 × diag(1,1,2) | raw voxel index -> code-contract um | 程序把 xyz 当像素坐标乘配置；不是文件自证的物理标定 | CONFIRMED_FROM_CODE [utils/rodent_vasculature/swc_io.py:185](</home/lzy/projects/ulm_3D_vascular/utils/rodent_vasculature/swc_io.py:185>) |
| 体素 spacing | [1, 1, 2] | um（配置） | 配置声明；代表 TIFF ResolutionUnit=1、无长度标定，实验依据 UNIT_UNVERIFIED | CONFIRMED_FROM_CONFIG [configs/swc_roi_generate.yaml:18](</home/lzy/projects/ulm_3D_vascular/configs/swc_roi_generate.yaml:18>) |
| 原始 SWC 半径 | 第六列原样保留 | 程序称 um；物理单位 UNIT_UNVERIFIED | 没有乘体素 spacing；不能因变量后缀认定原始单位 | CONFIRMED_FROM_CODE [utils/rodent_vasculature/swc_io.py:232](</home/lzy/projects/ulm_3D_vascular/utils/rodent_vasculature/swc_io.py:232>) |
| ROI 范围 | [80, 80, 120] | um（配置） | 采样几何参数，不是 PDMS 尺寸 | CONFIRMED_FROM_CONFIG [configs/swc_roi_generate.yaml:52](</home/lzy/projects/ulm_3D_vascular/configs/swc_roi_generate.yaml:52>) |
| 重采样间距 | 1.0 | um（配置） | 当前 smooth_centerlines=false；参数存在不等于每步启用 | CONFIRMED_FROM_CONFIG [configs/swc_roi_generate.yaml:25](</home/lzy/projects/ulm_3D_vascular/configs/swc_roi_generate.yaml:25>) |
| 重建半径系数 | 0.91 | dimensionless | 数值重建进料补偿，仅 H5 feed；原始半径保留 | CONFIRMED_FROM_CONFIG [configs/swc_stl_model_generate.yaml:21](</home/lzy/projects/ulm_3D_vascular/configs/swc_stl_model_generate.yaml:21>) |
| H5 第四列 | 2 * source_radius * 0.91 | diameter_um | 实际 H5 属性也写明 diameter_um | CONFIRMED_FROM_CODE [utils/cfd_lumen/ultraliser_backend.py:257](</home/lzy/projects/ulm_3D_vascular/utils/cfd_lumen/ultraliser_backend.py:257>) |
| Ultraliser 体素分辨率 | 6.0 | voxels/um | 重建工具参数，不是采集分辨率 | CONFIRMED_FROM_CONFIG [configs/swc_stl_model_generate.yaml:22](</home/lzy/projects/ulm_3D_vascular/configs/swc_stl_model_generate.yaml:22>) |
| STL 长度变换 | 1e-06 | um -> m | 只证明这条生成链的比例；不能推广到任意 STL | CONFIRMED_FROM_CODE [utils/cfd_lumen/ultraliser_qc.py:51](</home/lzy/projects/ulm_3D_vascular/utils/cfd_lumen/ultraliser_qc.py:51>) |
| 三角面目标边长 | 0.25913916380971913 | um | 派生表面重网格目标 | CONFIRMED_FROM_CONFIG [configs/cfd_surface_prepare.yaml:36](</home/lzy/projects/ulm_3D_vascular/configs/cfd_surface_prepare.yaml:36>) |
| 代表 ROI 半径范围 | [1.0, 2.7279] | um（保存元数据） | 派生直径为 2–5.4558 um；不是已验证的最小流体通道或 RBC 可通行性 | CONFIRMED_FROM_DATA [outputs/model_generate/ultraliser_anchor003274_20260907_192331/input/metadata.json:17](</home/lzy/projects/ulm_3D_vascular/outputs/model_generate/ultraliser_anchor003274_20260907_192331/input/metadata.json:17>) |
| 一维入口平均速度 | 0.7 | mm/s | 旧结构根节点边界假设；EXPERIMENTAL 未证明 | CONFIRMED_FROM_CONFIG [configs/cfd_preprocess.yaml:38](</home/lzy/projects/ulm_3D_vascular/configs/cfd_preprocess.yaml:38>) |
| 一维入口体积流量 | 7.693508475538942e-16 | m3/s | 历史 ROI 流量，不同于后期固定 CFD target | CONFIRMED_FROM_DATA [outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/qc/run_summary.json:24](</home/lzy/projects/ulm_3D_vascular/outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/qc/run_summary.json:24>) |
| 后期 CFD 体积流量 | 2.7369132390905703e-15 | m3/s | 旧 Tau1 baseline target；不是本轮测量 | CONFIRMED_FROM_CONFIG [configs/cfd_flow.yaml:65](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:65>) |
| 后期 CFD 质量流量 | 2.890180380479642e-12 | kg/s | 旧 CFD 参数 | CONFIRMED_FROM_CONFIG [configs/cfd_flow.yaml:64](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:64>) |
| 出口表压 | [14.544978101274268, 132.20454922317552, -13.700626673311461] | Pa | 延长段修正/数值基线；负表压不是绝对负压 | CONFIRMED_FROM_CONFIG [configs/cfd_flow.yaml:68](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:68>) |
| 结构叶压力参考 | 0.0 | Pa gauge | 结构叶零表压假设 | CONFIRMED_FROM_CONFIG [configs/cfd_preprocess.yaml:41](</home/lzy/projects/ulm_3D_vascular/configs/cfd_preprocess.yaml:41>) |
| 密度 | 1056.0 | kg/m3 | 旧均匀牛顿流体参数，不是当前样本测量 | CONFIRMED_FROM_CONFIG [configs/cfd_flow.yaml:55](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:55>) |
| 运动黏度 | 3.27e-06 | m2/s | 旧求解器/模型物性 | CONFIRMED_FROM_CONFIG [configs/cfd_flow.yaml:56](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:56>) |
| 动力黏度 | 0.00345312 | Pa s | 历史记录：rho * nu，不是新增 HemoCell 赋值 | CONFIRMED_FROM_DATA [outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/qc/run_summary.json:13](</home/lzy/projects/ulm_3D_vascular/outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/qc/run_summary.json:13>) |
| 体积黏度 | 2.18e-06 | m2/s | 旧 LBM 参数 | CONFIRMED_FROM_CONFIG [configs/cfd_flow.yaml:57](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:57>) |
| CFD 网格间距 | 2e-07 | m | 旧 Base 0.2 um 笛卡尔格，不是 HemoCell 网格决策 | CONFIRMED_FROM_CONFIG [configs/cfd_flow.yaml:46](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:46>) |
| CFD 时间步 | 2.038735983690112e-09 | s | 由 dx²/(6ν) 推导的旧 Tau1 数值步长 | INFERRED [utils/cfd_flow/validated_contract.py:60](</home/lzy/projects/ulm_3D_vascular/utils/cfd_flow/validated_contract.py:60>) |
| 数值压力偏置 | 3387510.7199999993 | Pa | rho * cs² * (dx/dt)²；不能解释为生理压力 | INFERRED [utils/cfd_flow/validated_contract.py:70](</home/lzy/projects/ulm_3D_vascular/utils/cfd_flow/validated_contract.py:70>) |
| 历史时间步 | 2.44140625e-08 | s | 明确 regression-only，不能与当前 Tau1 混用 | CONFIRMED_FROM_CODE [utils/cfd_flow/validated_contract.py:48](</home/lzy/projects/ulm_3D_vascular/utils/cfd_flow/validated_contract.py:48>) |
| NNE2 代表采集 XY 像素间距 | 1.16279069767442 | um/pixel | 显式采集 XML 元数据；只适用于该栈，不能用于 fMOST | CONFIRMED_FROM_DATA [vessel_model/T - NNE2/hana_stk/hana_stk/ZSeries-01212014-20xstack_site2-764/ZSeries-01212014-20xstack_site2-764.xml:31](</home/lzy/projects/ulm_3D_vascular/vessel_model/T - NNE2/hana_stk/hana_stk/ZSeries-01212014-20xstack_site2-764/ZSeries-01212014-20xstack_site2-764.xml:31>) |
| 通用 STL/NIfTI/Schmid/归档内部单位 | UNVERIFIED | UNIT_UNVERIFIED | 部分变量有 _um 后缀，但独立源单位/压力/时间证明不足 | UNVERIFIED [utils/schmid_pkl/loader.py:91](</home/lzy/projects/ulm_3D_vascular/utils/schmid_pkl/loader.py:91>) |


UNITS_STATUS：程序内部约定和部分文件单位可追溯；原始 fMOST 物理标定、原始半径、其他 STL/PKL/归档单位 `UNIT_UNVERIFIED`。既不全部宣布未知，也不把变量后缀当原始校准证明。

# 10. Coordinate system and transforms

| item | observed relationship | classification / evidence |
|---|---|---|
| origin | 活动 fMOST 链没有显式平移；局部 ROI 保留全局坐标。物理世界原点未知。 | CONFIRMED_FROM_CODE [utils/sampling/roi_extraction.py:89](</home/lzy/projects/ulm_3D_vascular/utils/sampling/roi_extraction.py:89>) |
| X/Y/Z | SWC 列 3/4/5 -> xyz；TIFF 数组为 zyx；解剖轴含义未证明。 | CONFIRMED_FROM_CODE [utils/rodent_vasculature/tiff_io.py:28](</home/lzy/projects/ulm_3D_vascular/utils/rodent_vasculature/tiff_io.py:28>) |
| axis direction / handedness | 缩放矩阵 diag(1,1,2) 不主动翻轴；输入是右手/左手、是否 LPS/RAS 均 COORDINATE_SYSTEM_UNVERIFIED。 | UNVERIFIED [utils/rodent_vasculature/swc_io.py:185](</home/lzy/projects/ulm_3D_vascular/utils/rodent_vasculature/swc_io.py:185>) |
| rotate / translate / center | 正式分区代码保留几何位置并写 translation_applied=False；可视化相机旋转不能当几何变换。 | CONFIRMED_FROM_CODE [utils/cfd_flow/geometry.py:179](</home/lzy/projects/ulm_3D_vascular/utils/cfd_flow/geometry.py:179>) |
| scale | fMOST 像素乘 spacing；后续 um -> m 乘 1e-6；这两次缩放不可漏做或重复做。 | CONFIRMED_FROM_CODE [utils/cfd_lumen/ultraliser_qc.py:51](</home/lzy/projects/ulm_3D_vascular/utils/cfd_lumen/ultraliser_qc.py:51>) |
| legacy image/world | main.py 的旧 STL/voxel 分支假设输入 LPS，write_nifti_mask 使用 X/Y 取反 affine 写 RAS 和 micron。输入 LPS 本身未被任意 STL 头证明。 | CONFIRMED_FROM_CODE [utils/io.py:38](</home/lzy/projects/ulm_3D_vascular/utils/io.py:38>) |
| NNE2 | TIFF stack 使用 vdb/XML spacing；XY 缩放经取整后保存 actual_x/actual_y；未推广为统一解剖坐标。 | CONFIRMED_FROM_CODE [utils/nne2/stack_io.py:120](</home/lzy/projects/ulm_3D_vascular/utils/nne2/stack_io.py:120>) |
| HemoCell lattice mapping | 轴选择、原点/边距、dx、STL 缩放策略、端口面到格点边界尚未建立。 | UNVERIFIED [configs/cfd_flow.yaml:46](</home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:46>) |


# 11. Boundary / inlet / outlet / wall labels

CONFIRMED_FROM_DATA：正式带标签 VTP 实读 73416 点、67262 cells，cell arrays 含 CellEntityIds、boundary_type_code、boundary_index、boundary_origin_code、boundary_origin、port_id、SurfaceRegionId/SurfaceRegion、RemeshEntityId。实体计数：wall ID1=67071；ID2=49；ID3=44；ID4=56；ID5=42。

boundary_manifest.csv 对应：inlet→ID4（CUT_PORT）；outlet_01→ID3（CUT_PORT）；outlet_02→ID5（CUT_PORT）；outlet_03→ID2（TRUE_TERMINAL）。四个端口 STL 实际存在。旧 CFD partition 从剩余实体推导 wall，并写五份米制 patch。STL 文件自身不能保存这些 VTP 属性，必须成套保留 CSV/JSON/VTP。

CONFIRMED_FROM_CODE / CONFIRMED_FROM_CONFIG：`port_transfer.py:50` 根据 local edge parent/child 方向产生 ASSUMED_INLET/OUTLET；TRUE_TERMINAL 被当作 ASSUMED_OUTLET；`measured_flow_direction=false`。因此这里既不是“只有无标签开口”，也不是“实验确认的入口”。HemoCell 如何使用这些面和 ID 属 UNVERIFIED，不能自动套用旧标签号。完整 CSV 记录见 `boundary_labels.json`。

# 12. Existing fluid / physics parameters

见 `physics_parameters.json`，每项有 old_value、old_unit、old_meaning、source_context、likely_origin、status。该表当前值全部是 `OLD_SOLVER_PARAMETER` 或其数值推导；没有将其标成 EXPERIMENTAL。

主要差别：一维 root=0.7 mm/s、ROI inlet Q=7.693508475538942e-16 m³/s；后期 Tau1 target=2.7369132390905703e-15 m³/s。二者属于不同历史尺度/边界契约，完整重标定来源待评审。rho=1056 kg/m³、nu=3.27e-6 m²/s、mu=0.00345312 Pa·s，不能直接作为 HemoCell plasma/RBC 模型参数。旧 P_ref≈3.3875 MPa 是数值压力偏置；历史 23622.320128 Pa 仅 regression reference。

# 13. Solver-independent components

CONFIRMED_FROM_CODE：SWC/ROI 数据结构、图连通性与 source-edge 映射、中心线/采样特征、网格拓扑、坐标变换、Ultraliser/VMTK 几何处理都没有依赖 HemoCell/Musubi 核心 API。独立于血流 solver 不代表无依赖或可无条件执行：VMTK/Ultraliser 有外部环境，某些 geometry loader 会在内存中合并点，repair/remesh 会生成改动几何。

INFERRED：优先 `utils/mesh/quality.py` 的无量纲拓扑函数；之后评审 `rodent_vasculature/swc_io.py`、`sampling/*`、`cfd_lumen/*`、`cfd_surface_prepare/*`。独立一维网络可以作为边界建模参考，但它是另一个物理模型，不能把结果视为实验真值。

# 14. Solver-specific components

CONFIRMED_FROM_CODE：`utils/cfd_flow/apes.py` 的 Lua 渲染与 WSL 启动、production/pipeline、restart_decode 的 LSB/PDF、Musubi one-step mass replay、D3Q19 link flux、固定 Tau1 数值参数，以及 patches/musubi 和 patches/seeder 均有旧 solver 耦合。

INFERRED：运行入口/输入适配未来应 REWRITE；重启、旧 core 补丁与固定参数 REFERENCE_ONLY；physical_port_flux、steady_state、steady_export 和 doc_visualize_v3 属 mixed，需 ADAPT 场命名/单位/边界契约。不能整套复制进 HemoCell，也不改 Palabos core。

# 15. Experimental / PDMS assets

见 `experimental_assets.json`。fMOST image/mask/SWC 和 NNE2 栈/XML 单列保护；NNE2 67 份采集 XML 有明确 micronsPerPixel。没有在已审阅文本、路径和代表元数据里确认 PDMS 几何、芯片尺寸、用户实测速度/流量/入口波形或 target-region 实验配准。`NOT_FOUND` 仅针对已声明搜索范围，不声称遍查 PDF 或 ZIP 内容。

Schmid 网络/RBC trajectories 文件有科研参考价值，来源和单位待确认，不能与采集图像或当前配置物性混为一类。

Schmid 数据族在清单中为 DATA_ONLY 候选，指后续可选择原始图结构；其压力/流量/RBC 轨迹等历史数值结果为 REFERENCE_ONLY。族级类别不授权整目录迁移。

# 16. Geometry-quality tools

| 检查/工具 | 实际位置与行为 | 迁移判断 |
|---|---|---|
| boundary/non-manifold/duplicate faces | utils/mesh/quality.py:97/108/138；VTK/NumPy，counts | ADAPT 最小只读接口 |
| degenerate / duplicate vertices | quality.py:80 与 cleanup.py:121 vtkCleanPolyData | ADAPT；必须分开检测和修复 |
| normals / connected components / bounds | cleanup.py:284 AutoOrientNormals；quality.py VTK connectivity、bounds/area/volume | ADAPT；面积和体积单位不得硬猜 |
| self-intersection | cfd_lumen/ultraliser_qc.py:174；trimesh 空间候选 + VTK triangle intersection；vmtk_qc/guarded_remesh 复用 | ADAPT；候选过滤和邻接排除有方法局限 |
| radius fidelity / minimum diameter | ultraliser_qc 的截面比较，input/metadata 的 source radius 统计 | ADAPT；不是完整最小流道或 RBC passage proof |
| repair / remeshing | pymeshfix cleanup；VMTK TPS、entity remesh、local_cut | ADAPT；未来单独授权几何变更 |
| voxelisation / skeleton | utils/voxel VTK image stencil、connectivity、skimage skeleton | ADAPT；不同图像轴/空间单位契约 |


本轮不执行这些处理或其测试。既有 QC 只作为 CONFIRMED_FROM_DATA 的历史结果，未重新证明网格质量。

# 17. External dependencies

`requirements.txt` 实际列出 NumPy、SciPy、h5py、VTK、PyVista、scikit-image、nibabel、Matplotlib、NetworkX、Jinja2、pymeshfix、psutil、pytest、Pillow、PyYAML、trimesh、Shapely、manifold3d、rtree；约束与 imports 见 `dependency_inventory.json`。静态分析不等于全部依赖已经可用。

CONFIRMED_FROM_DATA：本轮仅探测 `/usr/bin/python3` 的 module spec，numpy/h5py/PIL/yaml/vtk/scipy 可找到，pyvista/nibabel 不可找到；这不代表所有 Python 环境。没有安装依赖。VMTK Windows 环境路径在 WSL 不成立；APES 可执行路径仅 existence 检查，未运行。Ultraliser 是外部几何 C++ 工具，LBPM 是独立参考旧 solver，不能因二者在目录中而称它们都是活动主线。

CONFIRMED_FROM_CODE：Ultraliser 的 CMake 包含 OpenMP、TIFF、HDF5、Eigen3、GLM、FMT、ZLIB、BZip2；要求 CMake≥3.5、C++17，GNU 编译器实际版本检查为≥9.4（注释的 gcc 8.4 已过期）。这里只读构建声明，没有执行 CMake 或编译。

# 18. Migration candidates

`migration_inventory.json` 含 196 个模块/资产/配置条目，记录 source_path、purpose、entrypoint、solver_dependency、input/output、units、coordinate_system、dependencies、实验相关性、迁移类别、目标候选、证据和 unknowns；均为 DISCOVERED。

| class | 候选/范围 |
|---|---|
| KEEP_AS_IS | rodent_vasculature/geometry.py 的独立数学逻辑；sampling/feature_scaling.py 连同 ScalerState |
| ADAPT | SWC/ROI loader、mesh QC、几何/坐标、Ultraliser/VMTK 适配器、边界侧车和混合后处理 |
| REWRITE | APES 运行编排/输入渲染的功能接口，未来设计 HemoCell 应用层等价功能 |
| REFERENCE_ONLY | 旧 solver restart/PDF/core patches/Tau1 参数、历史验证及教学/文档结果 |
| DATA_ONLY | 经选择的原始图像/标注/ROI、单位明确的几何与标签侧车；原始大归档待评审 |
| DO_NOT_MIGRATE | Git 元数据、缓存、build、对象/二进制和打包运行时 |


目标侧 CONFIRMED_FROM_DATA：当前 HemoCell HEAD `5a410848bd5c57d5ae1c171112e78eab4a82e650`；已有 `M cmake/setup_googletest.cmake`、未跟踪验证 build/示例可执行文件及忽略的 build/tmp。它们全部是本轮开始前状态，不清理、不修复。根结构有 scripts/tools/examples/cases/data，没有 py_scripts/test_code。`cases/README.md` 自身称这些 cases 为潜在未验证/未完成案例。INFERRED 候选：只读几何工具可独立 tools/geometry_qc 或未来 py_scripts；真实求解应用才考虑 examples/cases；本轮不确定最终目录、不生成迁移代码。

# 19. DO_NOT_MIGRATE

CONFIRMED_FROM_DATA：.git（含 LFS 对象）、__pycache__/.pytest_cache/.ruff_cache/.codex_tmp/tmp、Ultraliser/build-wsl 的 CMakeFiles/CMakeCache/*.o/二进制、Windows BVLab-Annotation/_internal 运行时和 NNE2 安装器等不应进入新源码。source root 未发现 .venv/venv，但“数据目录”内部确实含打包库和 executable。

不得整目录丢弃 data/results/outputs：accepted Base restart、mesh/qval、物理指标、失败/未完成证据属于 REFERENCE_ONLY 科研凭证；仅区别“是否迁到新工程”，不是建议删除。零碎历史日志根据是否与验收绑定区分；本轮一律不清理。

# 20. Risks

| risk | level | reason |
|---|---|---|
| UNIT_RISK | HIGH | fMOST spacing 和原始半径单位主要来自代码/配置；泛用 STL 与代表 NIfTI 文件缺少单位证明；um/m 副本共存。 |
| COORDINATE_RISK | HIGH | 活动 xyz/zyx 与旧 LPS→RAS 分支不同；解剖轴、原点及 HemoCell 格点变换未确定。 |
| BOUNDARY_LABEL_RISK | HIGH | 标签存在，但角色是假设；STL 不携带 VTP 的 CellEntityIds/port_id，HemoCell 映射尚未实现。 |
| GEOMETRY_TOPOLOGY_RISK | MEDIUM | 历史 QC 通过不等于本轮重做拓扑验证，也不等于 HemoCell/RBC 分辨率、可通行性和边界可用性通过。 |
| SOLVER_COUPLING_RISK | HIGH | Musubi Lua/LSB、D3Q19 PDF、控制器、Tau1 固定契约不能直接用于 HemoCell；默认 FRESH_STEADY 实际被代码拒绝。 |
| EXPERIMENTAL_PROVENANCE_RISK | HIGH | 存在采集图像和 XML，但没有确认 PDMS、样本实测入口波形/流量、动静脉身份或原始 SWC 半径校准。 |
| DEPENDENCY_RISK | HIGH | VMTK Windows 路径在当前 WSL 不成立；requirements 与 /usr/bin/python3 环境不同；Git LFS 缺失；嵌套仓库 dirty。 |
| DATA_DUPLICATION_RISK | HIGH | 19 个 SWC 和 19 个 ROI 批次、不同模型/边界版本并存；同名原始 image 相同而 mask/SWC 不同，不能按名字去重。 |


# 21. Unknown / unverified items

- UNVERIFIED：fMOST 原始 SWC 半径是否已是物理 um，以及 [1,1,2] 间距的原始采集/论文对应证据。
- UNVERIFIED：各数据集的物理原点、轴方向、手性及其相互配准；不能把 legacy LPS/RAS 声明套到全部输入。
- UNVERIFIED：端口生理身份、实测入口流向/流量/波形、PDMS 通道及标定数据。
- UNVERIFIED：一维 ROI Q=7.6935e-16 与后期 target=2.7369e-15 的完整物理重标定依据；不能拼接不同阶段参数。
- UNVERIFIED：HemoCell 的几何单位契约、格点映射、边界实现与 RBC 所需有效流道分辨率。
- UNVERIFIED：多份兼容采样/模型输出中应冻结哪一个批次；null 自动解析不能代替选择决定。
- UNVERIFIED：VMTK 在本地 Linux 的可用环境及依赖版本组合；两个嵌套 Git 仓库 dirty 的具体源码/换行差异未深审。
- UNVERIFIED：通用源 STL/OBJ/NIfTI/Schmid PKL 和未解压 ZIP 的单位、授权、适用性及实验来源。
- UNVERIFIED：历史表面人工评审是否完成；本轮只读取现存 QC，不重跑自交/网格修复或 solver。

未完成事项的含义是“未来迁移的证据缺口”，不是已执行迁移的失败。本轮审计覆盖的文件/方法范围明确记录；最终只读一致性见 `READ_ONLY_VERIFICATION.json`。

# 22. Recommended first minimal migration unit

**FIRST_MINIMAL_MIGRATION_UNIT = READ_ONLY_TRIANGLE_TOPOLOGY_QC**（INFERRED；仅推荐）。

从 `utils/mesh/quality.py:32/97/108` 提取读取三角连接关系和边界边、非流形边、重复面计数的最小接口。输入是明确选择的三角 VTP 或 N×3 face IDs；输出是无量纲计数、输入身份和显式 unknown 单位/坐标元数据。它不改变 mesh，不运行 solver，不依赖实测流量，也不修改 HemoCell/Palabos core。

选择原因：当前代码和输入/输出已可静态确认，拓扑计数不使用物理坐标，因此不受尚未解决的单位/解剖轴影响；容易用四面体、缺面、非流形边、反序重复面验证。完整 `measure_mesh_quality` 的 _um 面积/体积字段应等待显式单位契约，暂不一并照搬。它也不能证明自交、法向、体积或 RBC 可通行性。

候选目标为独立 `tools/geometry_qc/` 或未来 `py_scripts/geometry_qc.py`，尚未选定或创建。详细接口/后续测试建议在 `first_minimal_migration_unit.json`；本轮没有执行这些测试或下一阶段。

NEXT_RECOMMENDED_ACTION = review audit results
