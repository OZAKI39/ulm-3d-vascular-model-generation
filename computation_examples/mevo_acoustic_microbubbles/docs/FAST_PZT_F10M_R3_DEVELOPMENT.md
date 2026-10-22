# F10M-R3：四PZT、双顶部孔与X两端封闭

本轮工程交付通过。正式四路装置声场、ON及导航均 **NOT_RUN**；未执行新的声波传播求解。已实际执行CPU几何/材料检查、四路kWaveArray源权重输入检查、GPU暴露壁面查询验证及16泡0.1 ms OFF完整M0短轨迹。顶部填充和自由液面、X端密封声学材料未知；training_ready=false、experiment_matched=false。

## 基线与独立执行

基于完整F10M提交a494cdddd71f0f7bbd09f4ef37e445a8b89d8548，保留F10L及F10M科学记录。独立分支dev/fast-pzt-f10m-r3-open-top-sealed-x-20261013，WSL工作区/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10m_r3_20261013/computation_examples/mevo_acoustic_microbubbles。远端/workspace/mevo_acoustic_microbubbles/f10m_r3_dual_port，使用既有隔离.venv，所有正式远端执行由Supervisor管理且来自.execution/<SHA>不可变源码。通过的GPU执行SHA为bf55520aabd715b138049e046fa0d65a8bd9a32d；绘图、CPU预检的执行SHA各自保留，不以最终文档提交替换旧运行身份。未更改旧求解器源文件或重跑旧Y管、B2及旧声场。

## 装置和实际连通性

PDMS为32.5×9×9 mm，原点位于几何中心；主通道沿X、ID1.6 mm、两侧PDMS3.7 mm，没有独立薄管壁。两个顶部孔ID1.5 mm，轴向Z，从Z=4.5 mm向下接入主通道，无底部贯穿和额外盲孔延伸。几何并集由同一implicit membership函数定义，材料标签、拓扑检查与显示使用同一来源；该函数不是精确SDF，不能直接减去泡半径作交汇处约束。

暂按端PZT发声面所在平面到相邻孔中心的X距离5.4 mm，得到孔中心(±10.85,0,4.5) mm，标PROVISIONAL_MEASUREMENT_REFERENCE，中心/边缘测距尚未确认。两孔镜像。50 µm几何拓扑网格六邻接洪泛得到1个连通空腔，主通道与两口属于同一分量；另检查轴向路径跨越交汇口，孔不穿底面。这个50 µm网格只是拓扑/显示用途，不作声学PPW或逐泡壁精度证据。

四枚SMD07T05R412WL有效陶瓷面直径7 mm：Y−中心(0,−4.5,0) mm法向+Y；Y+中心(0,+4.5,0)法向−Y；X−中心(−16.25,0,0)法向+X；X+中心(+16.25,0,0)法向−X；无Z源。两个端片都覆盖完整主通道端面。源面与孔体积的解析最短距离有实际可达的最近点见GEOMETRY.json；最小为相邻端面到孔体积的4.65 mm，未发现几何相交。检查对象是已知发声圆面和孔体积，未知PZT背衬/封装/电线未虚构建模。

## 顶部开放与X密封是不同边界

两X端机械边界永久封闭，球心须满足x∈[−16.25 mm+R,+16.25 mm−R]。闭端面作为零厚度机械接口记录，不为未知密封件创造实体厚度。材料、厚度、胶层、绝缘、超声透射/反射待确认，绝不自动设PDMS封头或理想刚性声学反射。原F10L开放端口实现保留以复现历史，但R3不继承StraightTubeSimulator，也不调用它的开放X端退出核。

顶部两口几何开放，既没有PDMS顶盖，也没有数值刚性顶盖。实际液位、孔内空气/液体、外部耦合介质及润湿条件未知。材料标签0外部未知、1PDMS文献参考、2主通道液体DEMO、3/4顶部空腔相态未知；boundary_id单独标两X机械密封面和两顶部开口。诊断材料网格100 µm，未知孔/外部/端密封声学节点保持NaN，现有spatial_medium适配器会拒绝，而非补水或透明端源。端面声学参数的未知状态优先于候选资源评估。

## 暴露壁面查询和GPU接入

新增coupling/pdms_cavity.py与fast_sim/pdms_cavity_engine.py，复用CylindricalLumen的坐标变换、GPUZeroFlow及FullParticleBatch完整M0。查询最近的暴露水平圆柱侧壁、两个竖孔侧壁、两X封闭圆盘及两组圆柱交线；被并集内部覆盖的圆柱片被排除。交线参数为x=c+a cosθ、y=a sinθ、z=√(R²−a²sin²θ)。CPU密采样1024个角度并逐个极小值括区优化，GPU128角度加40轮黄金分割，用于数值最近点而非伪称解析精确SDF。有限球体判断使用暴露边界距离减R，内部点/外部点符号来自并集成员关系。

在预测和每次PBD接触后调用新壁核；未改一次/二次Bjerknes力、原局部润滑、PBD、原0.1R墙修正及接触预算。修正过大时拒绝且不写回无效坐标。GPU 90个包含正负侧、交线内外和近壁位置的壁距与独立CPU参考最大差9.327e-17 m，有限球体投影最大差1.766e-14 m。另用1024个球面方向独立核对合法球体的几何成员关系。此一致性仅支持所测几何查询，不是局部微泡物理精度、时间步收敛或近壁水动力证明。

当前运动只准入z+R<1.5 mm的内部验证范围；1.5 mm是明确的数值范围界限，**不是实际液面位置**。孔附近验证假设该内部测试区域有液体，标INTERIOR_LIQUID_DEMO_ASSUMPTION_ONLY。超界返回131072并停止，绝不向下反射到虚构液面。真实自由液面运动、流动注入、残余流、声流和蒸发未实现。本轮对交汇口的改进是几何壁面查询已能使用，不是液面物理已验证。

## 实际检查、轨迹与成本

CPU 24项检查通过，预检实际1.166 s。四路CPU kWaveArray输入检查使用0.2 mm测试圆片、112³网格，不是7 mm真实装置缩尺预测；权重秩4、独立幅相/启停测试通过，BLI离边界19格，耗时0.777 s。没有调用声场run，也没有复用旧三路场。旧F10M的空间介质传播验证只作为已存在证据，不重算。

GPU 15项检查通过；16个固定ID、R=1.78 µm的确定性演示初态，覆盖主通道、两端附近、双交汇和孔内Z≤1.2 mm。真实OFF物理时长0.1 ms，dt1 µs，共100步、11帧；逐泡位置逐位不变，速度、背景梯度、声压及声力为零，重叠和墙修正为零。最小有限半径壁隙约1 nm，这是刻意布置的几何近壁数值检查，不是实验初态或物理近壁精度声明。合法半径/频率APECSS表SHA=5d53cb2ef42fa64591f45e3e7de5db106d6cbe74794c8a9abfa88a3354bec983；2.5 MHz标FREQUENCY_DEMO_2P5MHZ，零驱动为精确静止分支，不越过响应表下限外推。合法非零压力区间和有效叶节点保存在RESPONSE_SCOPE.json；没有匹配3 MHz表，ON拒绝。

通过的整组验证约7.912 s，含源输入、壁核编译及OFF；其中OFF积分/输出1.551 s，整卡采样峰值440 MiB。该耗时不含模块导入、先前失败、WSL绘图及归档，未作生产吞吐或长轨迹性能判断。此CPU oracle/GPU曲线搜索仍以正确性为本轮目标，未宣称它适合大规模生产。

正式四路声场和ON为NOT_RUN：液位/相态/外介质、端部密封声学结构、实际频率与电压声压标定尚不明确；候选3 MHz的64.8 GiB仍只是旧PDMS/水假设下的条件估计，若孔内空气或封装材料加入，最小声速、PPW、稳定步长和内存估计都需重新计算。没有对实际未知材料强行套用该网格。

## 失败保留和归档

CPU第一次预检的SLSQP线搜索失败，保存在failed_attempts/preflight_v1；改为轴对齐几何可证明下界且有可达最近点的解析距离，不放宽数值阈值。GPU validation_v1在Warp编译时发现动态循环布尔变量声明不符，已保存FAILURE和日志，改为bool(False)后validation_v2通过。两次源输入检查分别保留真实身份。首版图的文字重叠和不合适坐标自动偏移保存于figures/review_history_v1，修订只改变排版/坐标显示，没有改科学数据。

四组英文PNG/PDF已经实际打开检查，FIGURE_PROVENANCE、VISUAL_REVIEW、PDF_INSPECTION保存源文件与输出SHA。服务器重要文件31项/341548 bytes回传WSL并全部SHA一致；无未完成本轮求解。原F10L/F10M工作区保留。

## 文件与复现

- docs/PDMS_FOUR_PZT_OPEN_TOP_SEALED_ENDS.md：边界及接口说明。
- configs/straight_tube_exp1/f10m_r3/：明确几何、未知材料、内部测试范围、独立四路控制。
- results/straight_tube_exp1/f10m_r3/analysis/FINAL_GATE.json：工程通过和科学NOT_RUN分别记录。
- results/straight_tube_exp1/f10m_r3/validation_v2/：实际壁核比较、OFF轨迹、指标、响应表范围、GPU采样和源权重。
- figures/Fig_F10M_R3_1_apparatus、_2_cavity_connectivity、_3_source_port_geometry、_4_materials_and_admission：四组PNG/PDF。

在本项目目录用共享Python执行`PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.straight_tube_exp1.f10m_r3.verify`只读核对完整交付。重算必须另用结果目录和唯一Supervisor标签，提交源码后部署，不覆盖当前证据。

## 下一轮需要的实验信息

确认5.4 mm量到孔中心还是边缘、孔实际钻孔深度；两孔液位、液面形状/接触角、孔口上方和外部耦合介质；两X端实际密封材料、厚度、胶层、绝缘与PZT接触方式；PDMS品牌、配比口径、固化温度/时间、实验温度；通道液体；下一次单一频率、电压、四路相位及水听器声压传递标定；泡径和浓度。随后重新做材料/波长/资源准入、生成对应频率响应表和四个新底稿，再预登记短目标、对照和质量门槛。当前不开展正式导航优化。
