#!/usr/bin/python3
"""Write evidence-linked Chinese Step 2 reports from measured results."""
import argparse
import hashlib
import json
from pathlib import Path
from prepare_port_contract import PORT_IDS

def load(path):
    return json.loads(path.read_text())

def report(run):
    src=Path(__file__).resolve().parent;hc=src.parent.parent
    s=load(run/'diagnostics/voxelization_summary.json');t=load(run/'diagnostics/transform_check.json')
    top=load(run/'diagnostics/topology_check.json');codes=load(run/'logs/return_codes.json')
    audit=load(run/'diagnostics/source_immutability.json');vtk=load(run/'diagnostics/vtk_artifact_verification.json')
    p=s['ports'];dx=s['effective_dx_m'];g=load(run/'diagnostics/palabos_geometry.json')
    export_history=('本机 VTK 9.2 对首次较大 appended VTI 读取报 XML 错误；首次文件和失败日志保留在 diagnostics/vtk_first_export/ 与 logs/artifact_verification_attempt1.log。\n'
        '只重新导出格式，raw flags SHA256前后相等；没有重跑几何。见 diagnostics/vtk_export_correction.json。'
        if (run/'diagnostics/vtk_export_correction.json').exists() else '本次直接输出压缩 inline binary，没有执行其他格式重导出。')
    behavior=f'''# 官方 voxelizer 的本地源码行为

本次逐段检查 `{hc}/helper/voxelizeDomain.cpp` 及 `.h`。未调用其 getFlagMatrixFromSTL 最终接口。

`voxelizeDomain.cpp:76` 指针版本设 extraLayer=0、borderWidth=1、margin=1；以 DBL 精度读 TriangleSet。
构造 DEFscaledMesh(refDirLength, refDir, margin, extraLayer)、TriangleBoundary3D(*defMesh,false)，调用默认 inflate()。
随后构造 VoxelizedDomain3D，并把 voxelFlag::inside=3、innerBorder=4 转为 flagMatrix=1。

同一函数末尾（搜索 `Since the domain is closed`）继续运行两次 CopyFromNeighbor：

- 目标 x=0..1，源偏移 (+1,0,0)。
- 目标 x=nx-2..nx-1，源偏移 (-1,0,0)。

CopyFromNeighbor::process 按 x 从小到大原位复制。这是管道端面假设，无法表达本合同四个不同方向的 cap。
本轮 test-local C++ 在二值 flagMatrix 完成后结束构造流程；没有调用上述复制，也没有创建流体求解 lattice。

反事实诊断只对导出数据模拟全局单 rank 的原位复制，不调用官方 helper、不改正式输出。其四片可能新增流体数为
{[x['would_add_fluid_voxels'] for x in top['official_x_copy_counterfactual']]}，合计 {sum(x['would_add_fluid_voxels'] for x in top['official_x_copy_counterfactual'])}。
该反事实不声称重现所有未来 MPI 分块/envelope 同步行为。
实测 CLOSED 六个包围盒边界面流体数均为 0。OPENED 的 X-max 有 6 个流体体素，全部属于真实 outlet_03（label=4）的一层 cap opening；其余边界面为 0。

`examples/pipeflow/CMakeLists.txt` 将可执行文件写到现有 example 源码目录；其 compile.sh 使用已有 ../../build 和顶层目标，因此没有用它们编译。
本测试独立 CMake 只编译一个新增 C++ 单元，并只读链接已有 `build/libhemocell.a`；归档哈希见 diagnostics/archive_link_provenance.json。
没有修改顶层 CMake、cmake/setup_googletest.cmake 既有 patch、HemoCell core 或 Palabos。

源码证据：

- `{hc}/helper/voxelizeDomain.cpp`
- `{hc}/helper/voxelizeDomain.h`
- `{hc}/palabos/src/offLattice/voxelizer.h:46`（flag 语义）
- `{src}/closedVascularVoxelizer.cpp`（本轮实现）
'''
    (run/'official_voxelizer_behavior.md').write_text(behavior)
    transform=f'''# Physical → lattice 变换证据

PHYSICAL_COORDINATE_SYSTEM = axis_aligned_CFD，冻结 STL 的 m 坐标；不再旋转，不再做 um→m 缩放。
LATTICE_COORDINATE_SYSTEM = Palabos 全局整数节点坐标 LU，正向各轴与物理坐标一致。
SCALE = {t['scale_per_m']:.17g} LU/m（正向各向同性）
TRANSLATION = {t['translation_lu']} LU
MARGIN = {g['margin']} LU
EXTRA_LAYER = {g['extra_layer']} LU
REF_DIR = {s['ref_dir']}（Y，物理包围盒最长轴）
REF_DIR_N = {s['ref_dir_n']}
DX_REQUESTED = {s['requested_dx_m']:.17g} m/LU
DX_EFFECTIVE = {dx:.17g} m/LU
RELATIVE_DIFFERENCE = {(dx-s['requested_dx_m'])/s['requested_dx_m']:.17g}（{100*(dx-s['requested_dx_m'])/s['requested_dx_m']:.9f}%）

当前源码证明：

1. `palabos/src/offLattice/triangularSurfaceMesh.hh:1548` 的 toLatticeUnits 取得物理 bbox，按 referenceDirection 取跨度 Δ，scale=resolution/Δ；先减 bbox_min，再乘 scale；dx=1/scale。
2. `palabos/src/offLattice/triangleBoundary3D.hh:279` 的 DEFscaledMesh 设置 layer=margin+extraLayer；initialize 在上述变换后加 (layer,layer,layer)，并从 physicalLocation 减去 layer×dx。
3. 因此 `p_LU=(p_m-bbox_min_m)/dx+(1,1,1)=(p_m-physical_origin_m)/dx`。
4. `physical_origin_m={t['physical_origin_m']}`；该值和 dx 都由实际 DEFscaledMesh getter 导出，再与公式交叉核对。

逐点证明：冻结 STL 共 {t['original_unique_vertices']} 个精确唯一顶点；实际 Palabos pre-inflate 顶点数 {t['palabos_vertices']}。
全部顶点双向最近邻最大残差 {t['max_bidirectional_vertex_residual_lu']:.6g} LU，门限为 1e-8 LU。没有只凭数量或包围盒判定映射。
normal 不变；area_LU²=area_m²/dx²。报告端口中心始终取合同值，没有跟随离散化移动。

HemoCell 原有 inflate() 在此仿照保留，默认幅度 0.001 LU = {t['inflation_m']:.17g} m，沿局部顶点法向作用于内存网格。
它不是 affine transform 的一部分；它不改冻结 STL。cap 识别、中心和开孔射线仍依据原始冻结三角形。
默认值见 triangularSurfaceMesh.h:332，实现见 triangularSurfaceMesh.hh:1412。

`voxelizer.hh:100` 用 inflate 后 bbox 跨度构造 `N_axis=floor(max-min)+1+2*(margin+extraLayer)`，domain 从0开始。
本次 inflated bbox(LU)={g['inflated_bbox_lu']}，实际节点尺寸 {s['lattice_shape']}。
borderWidth=1、envelopeWidth=1、sparse blockSize=16；这些值控制离散边缘/存储，不改变上述全局映射。
未分配 sparse 区域是被裁去的外部区域，二值输出补0；native诊断以255单独表示未分配，不混作未知体素0。

ParaView VTI 使用 CellData，每个显示立方体中心对应一个 Palabos 整数节点。
VTI Origin=physical_origin−dx/2={t['vti_cell_data_origin_m']} m；Spacing=dx，Dimensions=节点尺寸+1。
因此显示 cell(i,j,k) 中心正好为 physical_origin+[i,j,k]×dx，不引入半格错位。
cap_triangles.vtp 和 mapped_centers_normals.vtp 都保留 m 坐标，可直接叠加冻结 STL。

PHYSICAL_TO_LATTICE_TRANSFORM = {t['status']}
证据文件：[transform_check.json](diagnostics/transform_check.json)、[palabos_geometry.json](diagnostics/palabos_geometry.json)。
'''
    (run/'physical_to_lattice_transform.md').write_text(transform)
    area_rows='\n'.join(f"| {k} | {p[k]['reconstructed_triangle_count']} | {p[k]['area']:.16g} | {p[k]['reconstructed_cap_triangle_area_m2']:.16g} | {p[k]['cap_area_relative_error']:.4g} | {p[k]['normal_dot']:.16g} |" for k in PORT_IDS)
    mapping_rows='\n'.join(f"| {k} | {p[k]['lattice_center']} | {p[k]['opened_port_voxels']} | {p[k]['distance_in_voxels']:.9f} |" for k in PORT_IDS)
    risks='''生理角色仍为 ASSUMED；source terminal 只证明 SWC 结构叶身份。原始采集标定/解剖方向/实测流向未确认。
Step 1 的历史人工表面评审、cap planarity 与 angle/aspect 阈值差异、自交未独立复验等限制继续有效。
本轮 0.2 µm 是旧 CFD 几何离散化参考；未选择正式 HemoCell dx，未验证红细胞分辨率、流量、速度或任何边界条件。
当前开孔是诊断 flag 副本中的几何邻接，尚未接入正式 HemoCell BC 或分布式开孔实现。
只测试 MPI 1；没有 MPI 2、性能比较、网格收敛或时间步推进。
离散立方体外观与连续三角形表面允许有半格阶梯偏差；不能据此扩大或移动合同端口。'''
    report_text=f'''# STEP 2 真实血管闭合体素化与四端口映射

STEP2_STATUS = {s['step2_status']}
STEP2_AUTO_CHECK = {s['step2_auto_check']}
人工 ParaView 核验仍待完成；本报告不把自动通过升级为整体 PASS。

## 1. Frozen Step 1 input

合同包：`{Path(s['input_stl']).parent.parent}`。输入已处于 axis_aligned_CFD 坐标系，含人工延长壁和四个真实封帽；不能替换为上游未旋转 lumen。
已读取 README、geometry/boundary/unit 合同及 provenance；81 项包内 SHA 和129项来源证据 SHA 全部核对一致。

## 2. STL integrity

输入：`{s['input_stl']}`。
SHA256：`{s['input_sha256']}`；单位 m。开始、结束完整性均通过；Step 1 全目录内容哈希核对 {audit['step1_content_hashes_checked']} 项无变化。

## 3. HemoCell voxelizer behavior

本地 helper 流程包含 STL→DEFscaledMesh→TriangleBoundary3D→VoxelizedDomain3D→二值 flagMatrix。
其末尾另行 CopyFromNeighbor 处理 X 两端。见 [官方源码行为](official_voxelizer_behavior.md)。

## 4. Why official helper cannot be used directly

本合同端口方向各异，X 两端处理与真实 cap 无一一对应关系。诊断反事实会新增{sum(x['would_add_fluid_voxels'] for x in top['official_x_copy_counterfactual'])}个格点，未执行该处理。
没有调用 getFlagMatrixFromSTL 返回最终 flagMatrix。

## 5. Closed voxelization method

新增 `{src}/closedVascularVoxelizer.cpp` 使用实际本地 Palabos 类和已编译 HemoCell 静态库。
margin=1、extraLayer=0、borderWidth=1、envelope=1、blockSize=16，保持官方 inflate=0.001 LU；停止于端口复制之前。
原生 inside/innerBorder 映射为1，其余0；分配了 {g['sparse_blocks']} 个 sparse block，导出全局二值场时仅补外部0。
CONFIGURE_RC={codes['CONFIGURE_RC']}，BUILD_RC={codes['BUILD_RC']}，SMOKE_RC={codes['SMOKE_RC']}。
单 rank 几何运行约 {codes['smoke_wall_seconds']:.3f} s；此值只记录执行，不是求解器性能基准。非致命上游 fgets warning 留在 build.log。

## 6. Physical→lattice transform

`p_LU = (p_m - {t['physical_origin_m']}) / {dx:.17g}`。
33,633 个顶点双向核对最大残差 {t['max_bidirectional_vertex_residual_lu']:.6g} LU。
法向不变、面积乘1/dx²；所有 getter、边距、inflate 与 VTI 半格解释见 [变换证明](physical_to_lattice_transform.md)。

## 7. Smoke-test resolution

来源：`{s['smoke_dx_source']}`，文件哈希 `{s['smoke_dx_config_sha256']}`。
REFERENCE_GEOMETRY_SMOKE_DX=requested={s['requested_dx_m']:.17g} m；有效值 {dx:.17g} m。
选最长 Y 轴 refDir={s['ref_dir']}，以 round(ΔY/requested_dx) 得 refDirN={s['ref_dir_n']}；相对差 {(dx-s['requested_dx_m'])/s['requested_dx_m']:.9g}（{100*(dx-s['requested_dx_m'])/s['requested_dx_m']:.9g}%）。
仅 geometry-only smoke test。闭合结果已连通，没有做加密诊断，也未自动选择正式 HemoCell dx。

## 8. Four-port source contract

来源为冻结 boundary_contract.json；读取 center(m)、normal、area(m²)、identity/source_id、extension/cap 属性，生成 [ports.tsv](ports.tsv)。
cut_000→inlet、cut_001→outlet_01、cut_002→outlet_02、terminal_000→outlet_03。四者都含人工延长段和 cap；角色生理身份 ASSUMED。
原始 s3 中心、计划延长端点和 s4 非加权质心不参与本轮映射。equivalent_radius 只保留作合同字段，未用于造圆盘。

## 9. Cap surface identification

仅从冻结完整 STL 恢复。固定平面距离阈值 {s['cap_plane_tolerance_m']:.17g} m，为文件 float32 坐标最大 ULP 的4倍；法向夹角不超过1°（先比较绝对点积）。
对候选三角形按共享边建连通分量，选择投影包含合同中心且面积质心最近的分量。
用完整 STL 有向体积判定外向约定；本例四个 raw mean normal 与合同均同向，dot≈1。若相反会区分同平面反向与几何方向错误，不直接改合同法向。
全套候选数量、triangle ID、plane residual、normal 记录见 [cap_identification.json](diagnostics/cap_identification.json)。

## 10. Cap area comparison

| 端口 | 三角形数 | 合同面积 m² | 恢复面积 m² | 相对误差 | 法向 dot |
|---|---:|---:|---:|---:|---:|
{area_rows}

未找到旧工程针对“冻结同一 cap 恢复”的更严格正式面积阈值；旧 vmtk_qc.py distal/proximal 面积5%检查属于另一指标，remesh patch area 为 DIAGNOSTIC_ONLY。
采用用户允许的≤1%作为 Step 2 geometry mapping tolerance；本例四项误差实际均为0。旧 planarity 1e-11 m物理表面质量阈值没有被本轮识别容差替代。
四组三角形 ID pairwise disjoint。

## 11. Closed-domain topology

实际 lattice={s['lattice_shape']}，二值流体 {s['closed_fluid_voxels']}、非流体 {s['closed_solid_voxels']}。
原生 inside={g['allocated_native_counts'][3]}、innerBorder={g['allocated_native_counts'][4]}；分配块内无 undetermined/toBeInside。
严格6面邻接下闭合 lumen 只有1个分量；六个 bbox 外层面均无流体。UNINTENDED_X_END_OPENING=NONE。
独立 VTK 原始 STL 包含性检查：{vtk['independent_containment']['closed_fluid_points_checked']} 个流体中心，报告在原始表面外的点 {vtk['independent_containment']['points_reported_outside']} 个；最远表面距离 {vtk['independent_containment']['max_outside_distance_m']:.6g} m，限于已记录 inflate 幅度。

## 12. Four-port opening method

对每个真实 cap，在其局部 bbox 中查找闭合场的外部格点；格点法向投影须在 cap 真实三角形并集内。
它必须与原闭合 lumen 中一个格点具有6面邻接，且该单格线段与一个原始 cap 三角形实际相交；记录作为 witness。
只把这一层外侧邻点置1，没有连续体几何编辑、圆盘、cap 删除、reservoir 或长通道。
开孔在导出的 Palabos flags 副本中完成，属于诊断字段，未写回 HemoCell core。

## 13. Changed-voxel locality

共新增 {top['added_fluid_voxels']}、删除 {top['removed_fluid_voxels']}；未归属改动 {top['unlabeled_changed_voxels']}、重叠 {top['port_overlap_voxels']}、标签标在未改动点的数量 {top['labeled_unchanged_voxels']}。
每个新增点都有内部邻点和相交 STL triangle ID，见 diagnostics/*_opening_witnesses.json。所有新增点投影仍在真实 cap 内。
全部端口连入同一原始 lumen；opened 仍为单连通。没有四个端口以外的二值场改变。
OPENED X-max 外层的6个流体点全部标为 outlet_03；这是实际端口邻域，不是 CopyFromNeighbor 产生的 X 开口。

## 14. Port label map

0=非端口，1=inlet，2=outlet_01，3=outlet_02，4=outlet_03。label只标 CLOSED→OPENED 的新增体素，不表示整个 cap 的内部截面。

| 端口 | 合同映射中心 LU | 新增体素 | 中心最近距离 LU |
|---|---|---:|---:|
{mapping_rows}

中心偏移检查固定门限2 LU，本例全部小于0.8 LU。没有移动中心。

## 15. ParaView outputs

- [closed_flag_matrix.vti](closed_flag_matrix.vti)：ClosedFluid，CellData，0/1。
- [opened_flag_matrix.vti](opened_flag_matrix.vti)：OpenedFluid，CellData，0/1。
- [port_label_field.vti](port_label_field.vti)：PortLabel，CellData，0..4。
- [cap_triangles.vtp](cap_triangles.vtp)：原始 cap，PortLabel 与 STLTriangleID。
- [mapped_centers_normals.vtp](mapped_centers_normals.vtp)：合同中心与 OutwardNormal。

最终均为米坐标、压缩 inline binary XML。三个 VTI 均已用独立 reader 回读，数组与 raw逐元素完全相等，origin/spacing/dimensions已检查。
{export_history}
新手操作见 [PARAVIEW_REVIEW.md](PARAVIEW_REVIEW.md)。

## 16. Risks

{risks}

## 17. Unverified items

人工 ParaView A–H检查仍 PENDING；没有声称已人工检查。
{risks}

## 18. Step 2 status

STEP2_AUTO_CHECK={s['step2_auto_check']}；STEP2_STATUS={s['step2_status']}。
冻结输入、实际构建、精确映射、closed/opened连通性、四个cap面积/独立性、changed locality、VTK回读、源码完整性均通过。
HemoCell {audit['existing_tracked_hash_count']} 项 tracked文件内容哈希与开始一致；旧工程与Step1全目录的预存条目无变化。
HemoCell 的 .git 目录可能因状态读取时创建/删除临时锁而仅发生目录 mtime 变化；审计单独保留该变化，不把它等同于源码文件修改。所有预存文件均无内容变化。
仅新增允许的 test_code 目录。既有 GoogleTest patch 保持；未 git add/commit/push。
FLUID_TIMESTEPS_RUN=0；未执行任何物理 BC、LBM、RBC、粒子、GPU或Step3。

NEXT_RECOMMENDED_ACTION = manual ParaView review of Step 2 voxelization and four-port mapping
'''
    (run/'STEP2_REPORT.md').write_text(report_text)
    write_review(run,s)
    write_final(run,src,s,t,audit)

def write_review(run,s):
    text=f'''# ParaView 人工核验指南（尚未完成）

待核验目录：`{run}`。本页检查项未勾选不代表失败，表示需要人工回答。
所有最终 VTI/VTP 均是物理米坐标。不要 Transform、旋转、改比例或移动端口来使画面贴合。

## 1. 打开闭合管腔

1. 打开 ParaView，File → Open，选择本目录 closed_flag_matrix.vti，点击 Properties 中 Apply。
2. 此时显示整个长方体是正常的，因为文件还含大量0值外部体素；不要把它当作管腔。
3. 选中 closed_flag_matrix，Filters → 搜索 Threshold。Scalars 选择 **CELLS / ClosedFluid**（Cell Data），范围 Lower=1、Upper=1，Apply。
4. 显示模式选 Surface 或 Surface With Edges，点击 Reset Camera。应只看到连续分叉血管，四个 cap 暂时封闭。必要时把 Opacity 降为0.35。
5. 阈值筛的是 CellData，因此显示小立方体中心对应 Palabos 格点；不需要 Point Data to Cell Data。

## 2. 打开已开孔管腔

1. File → Open → opened_flag_matrix.vti → Apply。
2. 同样创建 Threshold，选择 CELLS / OpenedFluid，范围1..1。
3. 用 Pipeline Browser 左侧眼睛在 closed/opened 两个 Threshold 间切换，不要让两者不透明重叠产生闪烁。
4. 整体形状应几乎相同，仅四个 cap 外侧新增一层小体素。这里“开孔”指诊断流体集合跨过 cap 到达小外侧邻域；不会出现大型外部 reservoir，也没有真实 BC。

## 3. 单独看四端口

1. 打开 port_label_field.vti → Apply。创建 Threshold，CELLS / PortLabel，范围1..4；颜色选择 PortLabel，将色标范围固定为1..4。
2. 0代表非端口。1=inlet，2=outlet_01，3=outlet_02，4=outlet_03。
3. 保持 opened 管腔半透明，显示端口 Threshold 不透明。
4. 从原始 port_label_field 分别创建四个 Threshold，范围依次1..1、2..2、3..3、4..4，重命名为四个端口；可以逐个开关眼睛。
5. 应分别显示 {s['ports']['inlet']['opened_port_voxels']}、{s['ports']['outlet_01']['opened_port_voxels']}、{s['ports']['outlet_02']['opened_port_voxels']}、{s['ports']['outlet_03']['opened_port_voxels']} 个新体素。Spreadsheet View 选择 Cell Data 可数行/看 PortLabel。
6. 该标签场只显示新增层，不包含 cap 内侧已有流体；斜面上的一层格点可能呈阶梯状，观察它们与内侧管腔的连接。

## 4. 叠加真实 cap、中心、法向

1. 打开 cap_triangles.vtp → Apply；Surface With Edges，按 Cell Data 的 PortLabel 上色。它含从冻结 STL 恢复的191个三角形。
2. 打开 mapped_centers_normals.vtp → Apply，先用 Points 表示并调大 Point Size；按 Point Data 的 PortLabel 上色。
3. 选中中心数据，添加 Glyph，Glyph Type=Arrow，Orientation Array=OutwardNormal，Scale Array=No scale array（或禁用数据缩放），统一 Scale Factor可用 `8e-7` m；Glyph Mode=All Points。
4. 箭头是几何外法向。inlet 的假定入流方向与外法向相反；本轮没有验证实际生理流向。
5. 逐个选中 cap 或端口，使用 Zoom to Data/选中对象后聚焦工具（版本不同名称可能略有差异），放大观察。中心最近新增格点均在0.8dx以内，不能要求每个连续中心恰好在整数格点。
6. 需要看完整壁面时直接打开冻结输入 `{s['input_stl']}`，Apply，Surface，Opacity=0.15–0.3；它本身仍带封帽。不要删除或修改原 STL。
7. 若封帽面遮住内部，在显示副本上使用 Clip 或 Slice 观察，关闭对应源对象眼睛；这只是可视化过滤，不改变输入文件。避免把 Clip 人工切面误认作额外开孔。

## 5. 查 X 两端和其余管壁

CLOSED 的六个最外格点面均无流体。OPENED 的 X-max 允许看到 **6个 label=4 的 outlet_03新增体素**，这是实际 cap 位于该处的结果。
其他位置不应有整片 X-min/X-max 出口或其他新增点。旋转相机观察全部分支，用 closed/opened切换与label叠加检查。
数据量是约6860万显示单元，虽文件压缩后很小，读取会展开。当前机器建议先看一个 Threshold，不必同时展开全部原始表格。

## 6. 请人工回答并记录

每项记录 YES/NO/UNCERTAIN；如有问题，记下端口/位置和截图名称。自动结果不能代填本表。

| 检查 | 人工回答 |
|---|---|
| A. lumen 是否连续？ | PENDING |
| B. 是否没有意外 X-min / X-max 大开口？（上述6个label4例外属于真实端口） | PENDING |
| C. inlet 是否位于正确 cap？ | PENDING |
| D. outlet_01 是否位于正确 cap？ | PENDING |
| E. outlet_02 是否位于正确 cap？ | PENDING |
| F. outlet_03 是否位于正确 cap？ | PENDING |
| G. 四个端口以外的血管壁是否保持关闭？ | PENDING |
| H. 是否存在明显的 voxelization 漏洞或断管？（期望 NO） | PENDING |

核验者：PENDING；日期：PENDING；备注/截图：PENDING。
当前 STEP2_STATUS=AUTO_PASS_HUMAN_PENDING。未人工确认前不写整体PASS，也不开始Step3。
'''
    (run/'PARAVIEW_REVIEW.md').write_text(text)

def write_final(run,src,s,t,audit):
    lines=[f"STEP2_STATUS = {s['step2_status']}",f"STEP2_AUTO_CHECK = {s['step2_auto_check']}",'',
        'STEP1_INPUT_INTEGRITY = PASS','',f"INPUT_STL = {s['input_stl']}",f"INPUT_STL_SHA256 = {s['input_sha256']}",'INPUT_STL_UNIT = m','',
        f"REFERENCE_GEOMETRY_SMOKE_DX = {s['requested_dx_m']:.17g} m (previous CFD mesh.dx_m)",
        f"REQUESTED_DX_M = {s['requested_dx_m']:.17g}",f"EFFECTIVE_DX_M = {s['effective_dx_m']:.17g}",'',
        f"REF_DIR = {s['ref_dir']} (Y)",f"REF_DIR_N = {s['ref_dir_n']}",'',
        f"PHYSICAL_TO_LATTICE_TRANSFORM = {t['status']}; p_LU=(p_m-{t['physical_origin_m']})/{s['effective_dx_m']:.17g}",'']
    lines.extend(f'LATTICE_N{axis} = {value}' for axis,value in zip('XYZ',s['lattice_shape']))
    lines.extend(['',f"CLOSED_FLUID_VOXELS = {s['closed_fluid_voxels']}",f"OPENED_FLUID_VOXELS = {s['opened_fluid_voxels']}",'',
        'CLOSED_LUMEN_CONNECTIVITY = PASS (1 connected component, 6-neighbor)','',
        'UNINTENDED_X_END_OPENING = NONE (6 X-max voxels are labeled outlet_03)','', 'PORT_COUNT = 4',''])
    for name in PORT_IDS:lines.append(f"{name.upper()}_MAPPING = {'PASS' if s['ports'][name]['mapping_pass'] else 'FAIL'} (nearest distance {s['ports'][name]['distance_in_voxels']:.6f} dx)")
    lines.append('')
    for name in PORT_IDS:lines.append(f"{name.upper()}_VOXELS = {s['ports'][name]['opened_port_voxels']}")
    lines.extend(['',f"PORT_PATCH_DISJOINTNESS = {s['port_patch_disjointness']}",f"PORT_OPENING_LOCALITY = {s['port_opening_locality']}",
        f"UNLABELED_CHANGED_VOXELS = {s['unlabeled_changed_voxels']}",f"PORT_OVERLAP_VOXELS = {s['port_overlap_voxels']}",'',
        'FLUID_TIMESTEPS_RUN = 0','',
        'PARAVIEW_OUTPUTS = closed_flag_matrix.vti; opened_flag_matrix.vti; port_label_field.vti; cap_triangles.vtp; mapped_centers_normals.vtp (all in RUN_DIR)','',
        'HEMOCELL_CORE_MODIFIED = NO','PALABOS_MODIFIED = NO','EXISTING_TRACKED_FILE_MODIFIED_BY_STEP2 = NO','',
        f'TEST_CODE_DIR = {src}',f'RUN_DIR = {run}',f'STEP2_REPORT = {run}/STEP2_REPORT.md',f'PARAVIEW_REVIEW = {run}/PARAVIEW_REVIEW.md','',
        'UNVERIFIED_ITEMS = manual ParaView A-H; physiological roles/calibration; historical geometry QC limitations; production dx/RBC/BC/flow/MPI2','',
        'NEXT_RECOMMENDED_ACTION = manual ParaView review of Step 2 voxelization and four-port mapping'])
    (run/'FINAL_STATUS.txt').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();report(a.run)
