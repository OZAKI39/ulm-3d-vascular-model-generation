# s3 实际调用链与冻结输入

入口完整读取：`/home/lzy/projects/ulm_3D_vascular/s3_cfd_1D_data_preprocess.py`；未执行任何源项目脚本。

| function | input | output | unit | boundary role | evidence |
|---|---|---|---|---|---|
| s3.main | 默认 configs/cfd_preprocess.yaml | CFDPreprocessConfig | um/Pa/m3/s | 调用入口 | [s3_cfd_1D_data_preprocess.py:45](</home/lzy/projects/ulm_3D_vascular/s3_cfd_1D_data_preprocess.py:45>) |
| load_cfd_preprocess_config | 严格 YAML 字段、PROJECT_ROOT 相对路径 | paths/selection/units/direction | 强制 um/Pa/m3/s | measured_flow_direction=false | [config.py:249](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/config.py:249>) |
| resolve_sampling_run / select_roi | 固定采样批次，anchor=3274 | 保存 ROIRecord；不重新采样 | 坐标/半径 um | 3 CUT_PORT、1 TRUE_TERMINAL | [pipeline.py:166](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/pipeline.py:166>) |
| resolve_rodent_run / load_matching_global_model | 当前 null 自动选最新兼容；实际冻结 input_manifest 记录 August run | 7419 节点/7418 边，验证 global_edges.csv | um | 原始 parent→child 结构 | [pipeline.py:176](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/pipeline.py:176>) |
| solve_global_flow | 原始半径、边、root=0.7 mm/s、叶0表压 | 全局 P/Q；不生成几何 | 内部 SI；输出 Pa/m3/s | LITERATURE_DERIVED_BASELINE_ASSUMPTION | [pipeline.py:185](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/pipeline.py:185>) |
| transfer_all_boundaries | ROI 边/全局节点映射与 P/Q | PortTransfer，保存角色/中心/法向 | um、无量纲、SI | cut 起点→ASSUMED_INLET；末点→ASSUMED_OUTLET；terminal需全局结构叶 | [port_transfer.py:245](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/port_transfer.py:245>) |
| build_boundary_geometry | 局部单入射边、中心、半径、角色 | 外法向、5D计划长度及端点 | um | 入口外法向与结构流向相反，出口同向 | [port_geometry.py:59](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/port_geometry.py:59>) |
| resolve_model_run / validate_model_run | 已生成模型和 QC，ROI/radius_scale/P95条件 | GeometryReference + 三个文件 SHA | um STL/VTP + m STL | 只作为表面参考，不改几何 | [pipeline.py:241](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/pipeline.py:241>) |
| write_boundary_package / write_extension_plan | 通过 readiness 的4个 PortTransfer | roi/boundary_conditions.json、CSV、port planes、extension plan | um、Pa、m3/s | 1 assumed inlet / 3 assumed outlets | [pipeline.py:289](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/pipeline.py:289>) |


实际代码先从全局图求一维解并传递边界，再校验已有 lumen 表面。STL 不是一维阻力求解输入；它是同 ROI 的只读几何参考。不能把流程简化成“读取 STL 后计算一维流量”。

本轮固定 s4 真正消费的 s3 包：`/home/lzy/projects/ulm_3D_vascular/outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628`。其 input_manifest 记录 sampling=`/home/lzy/projects/ulm_3D_vascular/outputs/sampling/20260825_133201_radius_plus_structure_k5`，rodent sample=`/home/lzy/projects/ulm_3D_vascular/outputs/rodent_vasculature/all_run_20260825_133152/samples/raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01`，model=`/home/lzy/projects/ulm_3D_vascular/outputs/model_generate/ultraliser_anchor003274_20260825_133350`。这些记录保留 E: Windows 前缀；对应 WSL 路径和逐文件 SHA 核对在 source_paths.json/hash_cross_checks.json，未修改原记录。

当前 s3 配置的 model_run/rodent_run 仍为 null；resolver 按目录 mtime/name 选最新兼容结果。当前文件中有 September 新批次，不能声称重跑当前 s3 一定产生被 s4 固定的 August 包。本轮以 s4 的显式输入及该包的历史 input_manifest 锁定这条链，未触发自动解析或重跑。

CUT_PORT 是 ROI 截边产生的边界；TRUE_TERMINAL 是保存原始 SWC 图的结构叶，不能升级为已证实的生理解剖终端。s3 的 role、原始中心/外法向/计划端点及源边 ID 完整保存在 boundary_contract.json。
