# RBC 几何适配硬门报告

RBC_GEOMETRIC_FIT = FAIL

本轮使用固定的原生 HemoCell RBC_FROM_SPHERE 网格。名义缩放 cbrt(50/90)=0.8220706914434901，名义直径 6.428592807 µm；实际顶点最大跨度 6.428992610 µm。网格为 642 顶点、1280 三角形，三个轴向跨度为 [6.428992605450912, 1.8856618638964933, 6.428992605450912] µm，最薄轴跨度约1.886 µm。未修改形状、网格分辨率或机械常数。

## 固定缩放的实际体积

原生离散网格实际体积=45.046330078320 µm³，面积=87.307010647669 µm²。与名义目标50 µm³的相对差异=-9.907339843%。50/90为用户给出的名义缩放关系，原生离散网格不是精确90 µm³的参考体积，因此不能把固定缩放后的几何声称为精确50 µm³。未为了消除该差异重新设定scale、半径或网格。本结果不是运行中的体积漂移。

## 真实网格放置搜索

闭合 STL 与 Step2 closed voxel mask 保持冻结。候选中心取既有lumen中每2 LU网格点，共22835个；131个确定性法向、每个0/30/60度旋转，共8974155个粗筛姿态。粗筛用真实网格的极值顶点与固定子集，再对候选检查全部642顶点；对48个候选按预先固定的8轮平移/旋转规则细化。它只调整预运行放置姿态，不调整细胞形状、尺寸或物理参数。

CASE A/C要求所有顶点有至少2dx的初始壁面间隙，并且RBC曲面与血管曲面没有相交。2dx是IBM核支撑对应的放置安全余量，不是允许穿墙的容差。穿墙容差为几何零：signed distance>0即越出闭合lumen。

没有找到可接受放置。最佳诊断候选仍有154/642个顶点在lumen外，最小间隙=-2.989411602672 LU=-0.597580930480 µm，检测到186个三角形接触对。这个候选仅供失败诊断，绝不是已冻结可运行的spawn。

此结论来自有界确定性搜索，不是对所有连续姿态的数学不可嵌入证明。当前没有获得可审核的安全初态，因此按本任务硬门停止；不能通过放宽间隙、缩小RBC或事后移动初始位置来宣称PASS。

## 独立壁面判据复核

以原始闭合STL的负signed distance为正向内壁间隙；闭合面包含端口cap，使用它对离开封闭lumen作保守筛查。保留native 0.001 LU inflate，不在诊断中扩大原始lumen。用既有inside格点及已知outside点验证sign；确认STL闭合且无非流形边。

另一个独立脚本以全曲面solid-angle winding number判定inside/outside，以独立nearest-triangle locator重算距离，以VTK mass properties重算体积/面积。与搜索代码一致，GEOMETRY_FINALIZER_IDENTITY=PASS。这只验证静态几何证据，不是SOLVER_FINALIZER_IDENTITY=PASS。

## ParaView 文件

geometry/FROZEN_LUMEN_LU.vtp 与 geometry/BEST_RIGID_PLACEMENT.vtp 使用同一LU坐标系，可共同打开查看。后者明确为失败诊断姿态；转换到物理坐标需使用contracts中的origin与dx。没有成功case，未伪造initial/mid/final运行帧。HUMAN_RBC_REVIEW=PENDING。
