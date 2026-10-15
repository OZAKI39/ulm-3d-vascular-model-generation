# 原生动力学实现与实际覆盖

只修改独立的四个BPM研究C++文件；pair_bpm_peri、verlet、fix_nve、fix_move原源码SHA不变。新内部style为LPS/VARIATIONAL/DYNAMIC all rho；要求SI、正密度、全组力换算、一个free组nve、仅允许move/property与该协调fix，禁止未审查的force modifiers。实际驱动明确不重叠的49/47/245边界与自由组。复制图上限4096节点不代表该规模已认证。

## 时间层

1. 原生free组NVE首次半kick和drift；边界move以连续物理time读取冻结解析曲线。
2. POST_INTEGRATE只在whichflag=1且非setup触发，以真实ntimestep为token；首次之后必须连续+1。old/new同坐标能差记D，只在新完全失效时更新owned键和special exclusions。
3. 此时早于Neighbor::decide。设置force_reneighbor及本步next_reneighbor；当前实测输入每步完整重建，pair在bond前计算但已看到新排斥拓扑。
4. 新损伤下pair/bond力密度汇总并reverse_comm，POST_FORCE只对owned乘同一V。setup也只转换本次新计算的初始力，不提交损伤。不改变原真实mass。
5. 原生NVE最后半kick。Python只读取接受结果/求积/停止，不做自由节点位置或速度推进。

## 状态与缓存

κ/d、原始键ID、固定m/r/taper、token、计数、D、接受坐标继续使用上一轮版本化bond restart。诊断求力次数不序列化，所以进程重启后它会重置；科学commit/token不重置。force/run0/output均不提交。静态显式C API在动态fix存在时返回−7。固定同坐标新状态缓存避免紧随commit的重复局部消元；d未变无需第二次消元。仅完整全图d=1可零能零力快速返回，部分损伤星邻域没有跳过被释放的内部自由度。

## 守恒与误差

native exported force为N，U与D为J；native epair乘V，native KE使用物理质量后已为J。实际记录全体K和K_free。夹持外功为−F_material的梯形路径功加ΔK_grips；外部冲量含ΔP_grips。储能与每步D由冻结离线模型在所有245保存态独立核对。已存终态3分量固定d能量差分最大相对误差1.3468e-9。

本分裂仅保证明确且一次提交；d改变时不宣称整体二阶准确。残差分解详见DYNAMIC_ENERGY_LEDGER.json和energy_residual_decomposition.npz。主工况Δd越界时停止，有限窗口dt/2改善不构成全程通过。

## 已实际执行与未覆盖

实际244个主物理步、一次同接受态restart续算、一次68步短dt/2、固定状态与只读试探。完整失效原键t=0；自然轨迹未闭合接触，未运行P4、2MPI动态迁移或生产规模。历史静态接触和2MPI证据不替代本轮未覆盖项。析构函数防护确保LAMMPS force先释放时不解引用悬挂material指针；正式测试前已修正。
