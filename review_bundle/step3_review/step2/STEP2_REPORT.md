# STEP 2 真实血管闭合体素化与四端口映射

STEP2_STATUS = PASS
STEP2_AUTO_CHECK = PASS
用户于 2026-09-13 在本会话确认：“人工 ParaView 检查通过”。结合既有自动检查 PASS，Step 2 总体状态更新为 PASS。详见 [人工确认记录](diagnostics/human_paraview_review.json)。

## 1. Frozen Step 1 input

合同包：`/home/lzy/projects/compre_output/step1/20260912_215759/geometry_contract`。输入已处于 axis_aligned_CFD 坐标系，含人工延长壁和四个真实封帽；不能替换为上游未旋转 lumen。
已读取 README、geometry/boundary/unit 合同及 provenance；81 项包内 SHA 和129项来源证据 SHA 全部核对一致。

## 2. STL integrity

输入：`/home/lzy/projects/compre_output/step1/20260912_215759/geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl`。
SHA256：`840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb`；单位 m。开始、结束完整性均通过；Step 1 全目录内容哈希核对 107 项无变化。

## 3. HemoCell voxelizer behavior

本地 helper 流程包含 STL→DEFscaledMesh→TriangleBoundary3D→VoxelizedDomain3D→二值 flagMatrix。
其末尾另行 CopyFromNeighbor 处理 X 两端。见 [官方源码行为](official_voxelizer_behavior.md)。

## 4. Why official helper cannot be used directly

本合同端口方向各异，X 两端处理与真实 cap 无一一对应关系。诊断反事实会新增205个格点，未执行该处理。
没有调用 getFlagMatrixFromSTL 返回最终 flagMatrix。

## 5. Closed voxelization method

新增 `/home/lzy/projects/hemocell_starter/test_code/step2_vascular_voxelization_smoke/closedVascularVoxelizer.cpp` 使用实际本地 Palabos 类和已编译 HemoCell 静态库。
margin=1、extraLayer=0、borderWidth=1、envelope=1、blockSize=16，保持官方 inflate=0.001 LU；停止于端口复制之前。
原生 inside/innerBorder 映射为1，其余0；分配了 305 个 sparse block，导出全局二值场时仅补外部0。
CONFIGURE_RC=0，BUILD_RC=0，SMOKE_RC=0。
单 rank 几何运行约 49.757 s；此值只记录执行，不是求解器性能基准。非致命上游 fgets warning 留在 build.log。

## 6. Physical→lattice transform

`p_LU = (p_m - [8.00903149195763e-05, 5.5701359566515553e-05, 0.00011143341591952393]) / 1.9989918081065344e-07`。
33,633 个顶点双向核对最大残差 1.30245e-13 LU。
法向不变、面积乘1/dx²；所有 getter、边距、inflate 与 VTI 半格解释见 [变换证明](physical_to_lattice_transform.md)。

## 7. Smoke-test resolution

来源：`/home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml:mesh.dx_m`，文件哈希 `1eecec545a5c021ab762a39dbe58dc2531a91083c223bdf4d23f76ea1b80a805`。
REFERENCE_GEOMETRY_SMOKE_DX=requested=1.9999999999999999e-07 m；有效值 1.9989918081065344e-07 m。
选最长 Y 轴 refDir=1，以 round(ΔY/requested_dx) 得 refDirN=494；相对差 -0.000504095947（-0.0504095947%）。
仅 geometry-only smoke test。闭合结果已连通，没有做加密诊断，也未自动选择正式 HemoCell dx。

## 8. Four-port source contract

来源为冻结 boundary_contract.json；读取 center(m)、normal、area(m²)、identity/source_id、extension/cap 属性，生成 [ports.tsv](ports.tsv)。
cut_000→inlet、cut_001→outlet_01、cut_002→outlet_02、terminal_000→outlet_03。四者都含人工延长段和 cap；角色生理身份 ASSUMED。
原始 s3 中心、计划延长端点和 s4 非加权质心不参与本轮映射。equivalent_radius 只保留作合同字段，未用于造圆盘。

## 9. Cap surface identification

仅从冻结完整 STL 恢复。固定平面距离阈值 5.8207660913467407e-11 m，为文件 float32 坐标最大 ULP 的4倍；法向夹角不超过1°（先比较绝对点积）。
对候选三角形按共享边建连通分量，选择投影包含合同中心且面积质心最近的分量。
用完整 STL 有向体积判定外向约定；本例四个 raw mean normal 与合同均同向，dot≈1。若相反会区分同平面反向与几何方向错误，不直接改合同法向。
全套候选数量、triangle ID、plane residual、normal 记录见 [cap_identification.json](diagnostics/cap_identification.json)。

## 10. Cap area comparison

| 端口 | 三角形数 | 合同面积 m² | 恢复面积 m² | 相对误差 | 法向 dot |
|---|---:|---:|---:|---:|---:|
| inlet | 56 | 7.819753007034155e-12 | 7.819753007034155e-12 | 0 | 1 |
| outlet_01 | 44 | 4.416140919927327e-12 | 4.416140919927327e-12 | 0 | 1 |
| outlet_02 | 42 | 4.315070751797983e-12 | 4.315070751797983e-12 | 0 | 1 |
| outlet_03 | 49 | 5.110621093985814e-12 | 5.110621093985814e-12 | 0 | 0.9999999999999999 |

未找到旧工程针对“冻结同一 cap 恢复”的更严格正式面积阈值；旧 vmtk_qc.py distal/proximal 面积5%检查属于另一指标，remesh patch area 为 DIAGNOSTIC_ONLY。
采用用户允许的≤1%作为 Step 2 geometry mapping tolerance；本例四项误差实际均为0。旧 planarity 1e-11 m物理表面质量阈值没有被本轮识别容差替代。
四组三角形 ID pairwise disjoint。

## 11. Closed-domain topology

实际 lattice=[493, 497, 280]，二值流体 182694、非流体 68423186。
原生 inside=113961、innerBorder=68733；分配块内无 undetermined/toBeInside。
严格6面邻接下闭合 lumen 只有1个分量；六个 bbox 外层面均无流体。UNINTENDED_X_END_OPENING=NONE。
独立 VTK 原始 STL 包含性检查：182694 个流体中心，报告在原始表面外的点 41 个；最远表面距离 1.98864e-10 m，限于已记录 inflate 幅度。

## 12. Four-port opening method

对每个真实 cap，在其局部 bbox 中查找闭合场的外部格点；格点法向投影须在 cap 真实三角形并集内。
它必须与原闭合 lumen 中一个格点具有6面邻接，且该单格线段与一个原始 cap 三角形实际相交；记录作为 witness。
只把这一层外侧邻点置1，没有连续体几何编辑、圆盘、cap 删除、reservoir 或长通道。
开孔在导出的 Palabos flags 副本中完成，属于诊断字段，未写回 HemoCell core。

## 13. Changed-voxel locality

共新增 497、删除 0；未归属改动 0、重叠 0、标签标在未改动点的数量 0。
每个新增点都有内部邻点和相交 STL triangle ID，见 diagnostics/*_opening_witnesses.json。所有新增点投影仍在真实 cap 内。
全部端口连入同一原始 lumen；opened 仍为单连通。没有四个端口以外的二值场改变。
OPENED X-max 外层的6个流体点全部标为 outlet_03；这是实际端口邻域，不是 CopyFromNeighbor 产生的 X 开口。

## 14. Port label map

0=非端口，1=inlet，2=outlet_01，3=outlet_02，4=outlet_03。label只标 CLOSED→OPENED 的新增体素，不表示整个 cap 的内部截面。

| 端口 | 合同映射中心 LU | 新增体素 | 中心最近距离 LU |
|---|---|---:|---:|
| inlet | [117.97642500628048, 9.269039858041094, 236.98437767350583] | 194 | 0.270522241 |
| outlet_01 | [172.02951952639518, 491.4707772181752, 273.8514889718654] | 85 | 0.550457701 |
| outlet_02 | [3.5505071178704655, 119.86441070663301, 6.528018106701496] | 93 | 0.737705531 |
| outlet_03 | [490.33398878144027, 413.64834078959814, 91.23076157835808] | 125 | 0.787709369 |

中心偏移检查固定门限2 LU，本例全部小于0.8 LU。没有移动中心。

## 15. ParaView outputs

- [closed_flag_matrix.vti](closed_flag_matrix.vti)：ClosedFluid，CellData，0/1。
- [opened_flag_matrix.vti](opened_flag_matrix.vti)：OpenedFluid，CellData，0/1。
- [port_label_field.vti](port_label_field.vti)：PortLabel，CellData，0..4。
- [cap_triangles.vtp](cap_triangles.vtp)：原始 cap，PortLabel 与 STLTriangleID。
- [mapped_centers_normals.vtp](mapped_centers_normals.vtp)：合同中心与 OutwardNormal。

最终均为米坐标、压缩 inline binary XML。三个 VTI 均已用独立 reader 回读，数组与 raw逐元素完全相等，origin/spacing/dimensions已检查。
本机 VTK 9.2 对首次较大 appended VTI 读取报 XML 错误；首次文件和失败日志保留在 diagnostics/vtk_first_export/ 与 logs/artifact_verification_attempt1.log。
只重新导出格式，raw flags SHA256前后相等；没有重跑几何。见 diagnostics/vtk_export_correction.json。
新手操作见 [PARAVIEW_REVIEW.md](PARAVIEW_REVIEW.md)。

## 16. Risks

生理角色仍为 ASSUMED；source terminal 只证明 SWC 结构叶身份。原始采集标定/解剖方向/实测流向未确认。
Step 1 的历史人工表面评审、cap planarity 与 angle/aspect 阈值差异、自交未独立复验等限制继续有效。
本轮 0.2 µm 是旧 CFD 几何离散化参考；未选择正式 HemoCell dx，未验证红细胞分辨率、流量、速度或任何边界条件。
当前开孔是诊断 flag 副本中的几何邻接，尚未接入正式 HemoCell BC 或分布式开孔实现。
只测试 MPI 1；没有 MPI 2、性能比较、网格收敛或时间步推进。
离散立方体外观与连续三角形表面允许有半格阶梯偏差；不能据此扩大或移动合同端口。

## 17. Unverified items

本次 ParaView A–H 检查已获用户整体通过确认；没有提供逐项独立答复、截图或实际检查时间，本记录不补造这些信息。以下物理和历史几何 QC 限制继续保留。
生理角色仍为 ASSUMED；source terminal 只证明 SWC 结构叶身份。原始采集标定/解剖方向/实测流向未确认。
Step 1 的历史人工表面评审、cap planarity 与 angle/aspect 阈值差异、自交未独立复验等限制继续有效。
本轮 0.2 µm 是旧 CFD 几何离散化参考；未选择正式 HemoCell dx，未验证红细胞分辨率、流量、速度或任何边界条件。
当前开孔是诊断 flag 副本中的几何邻接，尚未接入正式 HemoCell BC 或分布式开孔实现。
只测试 MPI 1；没有 MPI 2、性能比较、网格收敛或时间步推进。
离散立方体外观与连续三角形表面允许有半格阶梯偏差；不能据此扩大或移动合同端口。

## 18. Step 2 status

STEP2_AUTO_CHECK=PASS；STEP2_STATUS=PASS。
冻结输入、实际构建、精确映射、closed/opened连通性、四个cap面积/独立性、changed locality、VTK回读、源码完整性均通过。
HemoCell 847 项 tracked文件内容哈希与开始一致；旧工程与Step1全目录的预存条目无变化。
HemoCell 的 .git 目录可能因状态读取时创建/删除临时锁而仅发生目录 mtime 变化；审计单独保留该变化，不把它等同于源码文件修改。所有预存文件均无内容变化。
仅新增允许的 test_code 目录。既有 GoogleTest patch 保持；未 git add/commit/push。
FLUID_TIMESTEPS_RUN=0；未执行任何物理 BC、LBM、RBC、粒子、GPU或Step3。

NEXT_RECOMMENDED_ACTION = define Step 3 scope and physical boundary-condition contract
