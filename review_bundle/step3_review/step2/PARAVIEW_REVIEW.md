# ParaView 人工核验指南与通过记录

已核验目录：`/home/lzy/projects/compre_output/step2/20260912_225418`。用户于 2026-09-13 确认：“人工 ParaView 检查通过”。下表记录整体确认，未补造逐项回答或截图。
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
5. 应分别显示 194、85、93、125 个新体素。Spreadsheet View 选择 Cell Data 可数行/看 PortLabel。
6. 该标签场只显示新增层，不包含 cap 内侧已有流体；斜面上的一层格点可能呈阶梯状，观察它们与内侧管腔的连接。

## 4. 叠加真实 cap、中心、法向

1. 打开 cap_triangles.vtp → Apply；Surface With Edges，按 Cell Data 的 PortLabel 上色。它含从冻结 STL 恢复的191个三角形。
2. 打开 mapped_centers_normals.vtp → Apply，先用 Points 表示并调大 Point Size；按 Point Data 的 PortLabel 上色。
3. 选中中心数据，添加 Glyph，Glyph Type=Arrow，Orientation Array=OutwardNormal，Scale Array=No scale array（或禁用数据缩放），统一 Scale Factor可用 `8e-7` m；Glyph Mode=All Points。
4. 箭头是几何外法向。inlet 的假定入流方向与外法向相反；本轮没有验证实际生理流向。
5. 逐个选中 cap 或端口，使用 Zoom to Data/选中对象后聚焦工具（版本不同名称可能略有差异），放大观察。中心最近新增格点均在0.8dx以内，不能要求每个连续中心恰好在整数格点。
6. 需要看完整壁面时直接打开冻结输入 `/home/lzy/projects/compre_output/step1/20260912_215759/geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl`，Apply，Surface，Opacity=0.15–0.3；它本身仍带封帽。不要删除或修改原 STL。
7. 若封帽面遮住内部，在显示副本上使用 Clip 或 Slice 观察，关闭对应源对象眼睛；这只是可视化过滤，不改变输入文件。避免把 Clip 人工切面误认作额外开孔。

## 5. 查 X 两端和其余管壁

CLOSED 的六个最外格点面均无流体。OPENED 的 X-max 允许看到 **6个 label=4 的 outlet_03新增体素**，这是实际 cap 位于该处的结果。
其他位置不应有整片 X-min/X-max 出口或其他新增点。旋转相机观察全部分支，用 closed/opened切换与label叠加检查。
数据量是约6860万显示单元，虽文件压缩后很小，读取会展开。当前机器建议先看一个 Threshold，不必同时展开全部原始表格。

## 6. 人工确认记录

确认来源为用户消息，不是自动检查推定。用户给出整体通过确认；下表各项纳入本次确认范围，未提供独立 YES/NO 答复。

| 检查 | 人工回答 |
|---|---|
| A. lumen 是否连续？ | PASS（用户整体确认） |
| B. 是否没有意外 X-min / X-max 大开口？（上述6个label4例外属于真实端口） | PASS（用户整体确认） |
| C. inlet 是否位于正确 cap？ | PASS（用户整体确认） |
| D. outlet_01 是否位于正确 cap？ | PASS（用户整体确认） |
| E. outlet_02 是否位于正确 cap？ | PASS（用户整体确认） |
| F. outlet_03 是否位于正确 cap？ | PASS（用户整体确认） |
| G. 四个端口以外的血管壁是否保持关闭？ | PASS（用户整体确认） |
| H. 是否存在明显的 voxelization 漏洞或断管？（期望 NO） | PASS（用户整体确认） |

核验者：用户（未提供姓名）；确认日期：2026-09-13；实际检查时间及截图：未提供。

原始确认文字、记录时间与核验产物哈希见 [human_paraview_review.json](diagnostics/human_paraview_review.json)。
当前 STEP2_STATUS=PASS；STEP2_AUTO_CHECK=PASS；HUMAN_PARAVIEW_REVIEW=PASS。此更新仅记录 Step 2 验收，未执行 Step 3。
