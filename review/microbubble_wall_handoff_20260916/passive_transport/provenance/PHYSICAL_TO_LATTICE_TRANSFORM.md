# Physical → lattice 变换证据

PHYSICAL_COORDINATE_SYSTEM = axis_aligned_CFD，冻结 STL 的 m 坐标；不再旋转，不再做 um→m 缩放。
LATTICE_COORDINATE_SYSTEM = Palabos 全局整数节点坐标 LU，正向各轴与物理坐标一致。
SCALE = 5002521.7509380905 LU/m（正向各向同性）
TRANSLATION = [-400.6535424246619, -278.64726278831756, -557.4480869187494] LU
MARGIN = 1 LU
EXTRA_LAYER = 0 LU
REF_DIR = 1（Y，物理包围盒最长轴）
REF_DIR_N = 494
DX_REQUESTED = 1.9999999999999999e-07 m/LU
DX_EFFECTIVE = 1.9989918081065344e-07 m/LU
RELATIVE_DIFFERENCE = -0.00050409594673276475（-0.050409595%）

当前源码证明：

1. `palabos/src/offLattice/triangularSurfaceMesh.hh:1548` 的 toLatticeUnits 取得物理 bbox，按 referenceDirection 取跨度 Δ，scale=resolution/Δ；先减 bbox_min，再乘 scale；dx=1/scale。
2. `palabos/src/offLattice/triangleBoundary3D.hh:279` 的 DEFscaledMesh 设置 layer=margin+extraLayer；initialize 在上述变换后加 (layer,layer,layer)，并从 physicalLocation 减去 layer×dx。
3. 因此 `p_LU=(p_m-bbox_min_m)/dx+(1,1,1)=(p_m-physical_origin_m)/dx`。
4. `physical_origin_m=[8.00903149195763e-05, 5.5701359566515553e-05, 0.00011143341591952393]`；该值和 dx 都由实际 DEFscaledMesh getter 导出，再与公式交叉核对。

逐点证明：冻结 STL 共 33633 个精确唯一顶点；实际 Palabos pre-inflate 顶点数 33633。
全部顶点双向最近邻最大残差 1.30245e-13 LU，门限为 1e-8 LU。没有只凭数量或包围盒判定映射。
normal 不变；area_LU²=area_m²/dx²。报告端口中心始终取合同值，没有跟随离散化移动。

HemoCell 原有 inflate() 在此仿照保留，默认幅度 0.001 LU = 1.9989918081065343e-10 m，沿局部顶点法向作用于内存网格。
它不是 affine transform 的一部分；它不改冻结 STL。cap 识别、中心和开孔射线仍依据原始冻结三角形。
默认值见 triangularSurfaceMesh.h:332，实现见 triangularSurfaceMesh.hh:1412。

`voxelizer.hh:100` 用 inflate 后 bbox 跨度构造 `N_axis=floor(max-min)+1+2*(margin+extraLayer)`，domain 从0开始。
本次 inflated bbox(LU)=[[0.999129717568173, 0.9990006572477488, 0.99903105551689], [491.17320385770046, 495.0008923850483, 278.6349897831459]]，实际节点尺寸 [493, 497, 280]。
borderWidth=1、envelopeWidth=1、sparse blockSize=16；这些值控制离散边缘/存储，不改变上述全局映射。
未分配 sparse 区域是被裁去的外部区域，二值输出补0；native诊断以255单独表示未分配，不混作未知体素0。

ParaView VTI 使用 CellData，每个显示立方体中心对应一个 Palabos 整数节点。
VTI Origin=physical_origin−dx/2=[7.999036532917097e-05, 5.560140997611023e-05, 0.0001113334663291186] m；Spacing=dx，Dimensions=节点尺寸+1。
因此显示 cell(i,j,k) 中心正好为 physical_origin+[i,j,k]×dx，不引入半格错位。
cap_triangles.vtp 和 mapped_centers_normals.vtp 都保留 m 坐标，可直接叠加冻结 STL。

PHYSICAL_TO_LATTICE_TRANSFORM = PASS
证据文件：[transform_check.json](diagnostics/transform_check.json)、[palabos_geometry.json](diagnostics/palabos_geometry.json)。
