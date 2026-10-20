# CLOT_SPATIALLY_HETEROGENEOUS_VARIATIONAL_LPS_V1

运行：`20261017T214917Z`。标签 `SPATIALLY_HETEROGENEOUS_EXPLORATORY_UNCALIBRATED`。

## 本轮结论

在本次单一冻结材料场与同一开环剪切驱动下，H未达到至少5个20–200节点凝聚性游离块的完整探索目标。必须区分损伤路径改变、失效键数量和真实中等块体形成，不能用前两者代替后者。

U：末态20–200节点块 0；历史最多 0；最大自由块 2 节点；细屑/脱落质量 100.0000%；实际净外功 13.940604 nJ。；H：末态20–200节点块 0；历史最多 0；最大自由块 3 节点；细屑/脱落质量 95.9459%；实际净外功 13.048879 nJ。

## 冻结模型、场与能量实现

1792材料点、336固定锚点、95666原始键、2736原始预切；参考X、m、r0、taper、质量及邻域不变。dx=50 µm、horizon=150.5 µm、G=940 Pa、K=100000/3 Pa、ρ=1000 kg/m³。bulk s0=0.15；interface s0/sc=0.140625/0.16875；固定支架卡不变。所有材料属性以原始atom/bond ID对应。需要明确的初始化细节：未预切键初始κ=0、d=0；原始预切d=1保持。现有初始化用κ=sc作为预切完全失效标记，因此H预切键的κ标记随其局部sc变化，数值上不与U逐键相同。这不是新增变形或新增预切，已完全失效键仍无材料力；但不能声称两组所有初始κ数组逐值完全相同。

GSTools 1.7.0；Gaussian z方差1，各向同性，显式rescale=1，C_z(r)=exp[-(r/300 µm)²]；RandMeth 4096模式，seed=20261017。q=tanh(z−c)，c使材料参考体积加权q均值为0。实际q范围[-0.954911152,0.925324541]，标准差0.507946320。单一冻结实现，没有结果后换种子。

对全部原始bulk键使用qij=(qi+qj)/2、Δsij=(0.18−0.15)exp(0.30qij)/Z、scij=0.15+Δsij。Z=1.01884313127；参考权重w=h*r0*V²*(1/mi+1/mj)，与完好参考水静pair刚度乘r0²成比例，包括原始bulk预切。sc范围0.172129358950–0.188831493601，加权平均0.180000000000；排除预切后的均值0.179997716424。归一化软化区间不等于归一化所有开裂路径耗能。

完整能量仍为每个原始端点star的正定局部弹性/体积能与串联柔度罚能之和，对内部弹性伸长精确消元。固定d时，U(x,d)、−∂U/∂x及力装配公式未改；sc只进入原生不可逆演化 d(κ)=0 (κ≤s0)、sc(κ−s0)/[κ(sc−s0)] (s0<κ<sc)、1 (κ≥sc)。κ=max历史正伸长。D按固定当前x的U(x,d_old)−U(x,d_new)累加，数值R单列，不混入D。没有在节点力上乘随机因子。κ驱动本构没有显式速率变量；本轮惯性/拖曳响应不能视为Liu实验内禀断裂速率依赖的复现。

新增complete_by_bond与原始bond ID冻结接口、restart V4；读旧V2/V3时回退typed sc。原有onset字段保留。设置接口必须在初始history导入/commit前调用，拒绝缺失/重复ID、非法sc及界面属性改动。输入逆序与atom排序后仍保持对应。隔离库及其工厂、协调器依赖重编译；旧库保持原哈希。

新库SHA `f2604a7f1023c75e6d527af463f93eb030b12cfb97fd9d2b0ffa8eb7110a4adb`；冻结SHA `a32edb12faaa1a57997b4e6dac3873081758977a40084f766f2df9fb48795d28`。原surface实现逐字节复用；积分仍为原生VV、速度相关末半步隐式回调、BASE=0.25 µs、最小预测子步0.03125 µs、Δd≤0.1及能量门槛1%。预测拒绝发生于native commit之前；提交后异常会保留failed_unaccepted而非声称已回滚。

## 实际材料场结构与尺度

线性梯度R²=0.122504；26邻接局部极值统计{'maxima': {'all_nodes': 4, 'interior_nodes': 1}, 'minima': {'all_nodes': 9, 'interior_nodes': 0}}；q>0与q<0分别1、1个连通区域。试样材料中心跨度650×350×750 µm；lc/dx=6，lc/horizon=1.993355。正负区域连通不等于只有两个组织，也不证明各局部极值统计独立。没有测得足够多的独立材料块；有限域和非局部作用可能削弱局部属性对分离路径的辨识。

## 规定路径耗能（离线，不是native动力学，也不是Gc）

完整原图在参考z=0两侧作±opening/2刚体平移；没有新增初始预切。0–180 µm，共721规定状态。初始刚度相同；第一非零损伤网格点7.75 µm并非精确起裂点。该点已经进入软化，反力不要求相同；严格起裂前的弹性与统一κ=s0阈值保持。

|离线路径材料卡|完整分离损伤耗能 nJ|末端功积分余额 nJ|
|---|---:|---:|
|uniform|7.83823688|-0.09501564|
|spatial_actual|8.05785554|-0.09644179|
|minimum|7.26202485|-0.09122786|
|median|7.90161442|-0.09545048|
|maximum|8.51807596|-0.09937412|

3552条未预切跨面原始键最终全部失效；损伤能量增量非负。耗能仍有有限规定步长误差：外功梯形积分与储能+耗能约1.2%差异，不能将上述值称精确连续极限或GC_EFFECTIVE。最小/中位/最大sc三张卡仅为离线诊断，不是新增正式剪切工况。

## 同一开环加载与实际成本

从上一轮checkerboard接受记录提取14840个幅值结点，含t=0时3200 Pa，按时间分段线性插值；单独乘原2 ms升、2.4 ms降、2.8 ms停驱包络。与旧FP64回调幅值最大误差7.97e-10 Pa。两个正式case共享同一只读文件，完全取消等功反馈。250 µm棋盘空间模板、当前位置/法向、被动拖曳与开放裂面代理不变。实际两组功允许不同。这里是数千Pa机械需求剪切代理，不是已验证的真实微泡声流。

|工况|净外功 nJ|驱动功 nJ|拖曳功 nJ|实际峰值剪切 Pa|最大能量残差 %|wall s|CPU s|峰值RSS MiB|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|uniform|13.940604|14.368617|-0.428013|5901.122996|0.296424|2166.185|2165.470|334.406|
|heterogeneous|13.048879|13.437539|-0.388660|5789.717824|0.286677|2139.171|2138.554|336.707|

U：最大面片速度0.974892 m/s；规定外缘速度代理峰值147.764005 m/s；Re_h峰值3694.100；初始受力面积1.680000 mm²，按首次受力面积累计的独特面片面积2.395953 mm²；H：最大面片速度0.884124 m/s；规定外缘速度代理峰值145.439359 m/s；Re_h峰值3635.984；初始受力面积1.680000 mm²，按首次受力面积累计的独特面片面积2.331723 mm²。数千量级Re_h不支持把该薄层滑移代理解释为低Re Stokes预测；规定外缘速度不是CFD求解或实验测得的微泡声流。各时刻有效受力面积保存在LOADING_HISTORY.csv。

完整账本：W_drive+W_drag=ΔU+ΔK+D+ΔUc−Jc+R。实际接受端点力平均点乘原生位移；Jc单列新接触排除解除产生的能量跃变。残差分母是历史最大能量尺度，不是简单末态W。所有接受行、每次损伤变化与真实分量成员保留；完整FP64坐标/速度/κ/d每25 µs及事件保存，未声称每接受步均保存完整x/v。

## 必须回答的28项科研问题


1. 已实际接入原生不可逆更新及完整变分能量路径。sc参与κ→d；力和U通过实际d变化响应，未改最终力乘数。


2. 是。初态参考文件逐字节沿用；原始节点、95666键、2736预切及m不变；本轮检查和正式末态原生图核验通过。


3. 不是简单线性梯度（R²=0.1225），有平滑三维变化；极值和正负连通域如上，经验相关函数在Fig2。


4. 尚不能证明足够多统计独立材料区。大多数极值位于边界，正负场各连通，颜色数量不是组织区域数量。


5. bulk sc最小0.172129358950、最大0.188831493601、参考加权均值0.180000000000。


6. 保持。q=0新旧库三种状态逐位一致；H和U完好能量/弹性一致，bulk s0全部0.15，界面不改。


7. 保持固定历史下的力—能量一致性；三个代表状态最大相对有限差分误差2.702e-08，且异质低/高q节点均覆盖。力矩、总力及独立完整能量核对通过。


8. 保持原生不可逆规则。未预切键初始κ/d一致；预切d=1相同，κ=sc初始化标记的数值差异如上披露。卸载不恢复、κ不下降、未提交试探不改历史；排序改变2126个atom槽位后ID/参数保持，restart V4保留sc、κ/d和D。未重跑完整MPI矩阵。


9. 是，同一冻结开环幅值、空间函数、包络和表面算法。不同当前表面和速度导致实际牵引及功不同，属于材料响应差异。


10. U：净外功13.940604 nJ，驱动14.368617，拖曳-0.428013 nJ；H：净外功13.048879 nJ，驱动13.437539，拖曳-0.388660 nJ


11. U：0.796500 ms，原始 bond ID [32640]，参考中点 [[-0.000125, 0.00015000000000000001, 0.000375]] m；H：0.796500 ms，原始 bond ID [32640]，参考中点 [[-0.000125, 0.00015000000000000001, 0.000375]] m。首次原始键相同=True。


12. U：0.838250 ms，原始 bond ID [32640]，参考中点 [[-0.000125, 0.00015000000000000001, 0.000375]] m；H：0.834750 ms，原始 bond ID [60594]，参考中点 [[0.000125, -0.00015000000000000001, 0.000375]] m。首次原始键相同=False。


13. U：失效键共享端点网络出生56次，最多同时11区；H：失效键共享端点网络出生58次，最多同时14区。这是损伤区域代理，不能自动解释为独立裂纹起源。


14. 真实三维连续裂面贯通尚未由完整连续表面证明。自由连通分量证明其与锚/主体之间的所有原始存活材料路径已切断，但小点链分量不能当作中等三维裂块。辅助开放裂面图和failed-edge区数均不充分。


15. U：自然bulk失效15010，interface失效0；H：自然bulk失效13919，interface失效0


16. U：末态20–200块0，历史最多0，目标质量/初始0.0000%；H：末态20–200块0，历史最多0，目标质量/初始0.0000%


17. U：末态最大2节点，历史最大3，最大/脱落质量2.5000%；H：末态最大3节点，历史最大4，最大/脱落质量4.0541%


18. U：1–2节点细屑/脱落质量100.0000%，总脱落/初始4.4643%；H：1–2节点细屑/脱落质量95.9459%，总脱落/初始4.1295%


19. uniform：全部二次分裂3次，母块≥20节点0次；heterogeneous：全部二次分裂5次，母块≥20节点0次。细屑分裂不称凝聚性中等块二次断裂；内部存活键、birth、寿命和体积等效直径保存在真实拓扑历史。


20. uniform：q与首次软化时间的条件Spearman=0.35549488807305685，q与末态d=-0.0428782806674881；heterogeneous：q与首次软化时间的条件Spearman=0.31303279182584487，q与末态d=-0.02638373820382995。低/高q四分位以全部合格原始bulk键为分母；U中q只是参考协变量。


21. 只能描述关联。只有一个随机实现和一对对照；位置、原始切口、暴露表面、加载及非局部耦合均会混杂。没有独立样本数、显著性或组织组成因果结论。


22. 存在：lc/δ=1.9934，lc/dx=6。试样Y方向仅350 µm，约1.17lc；非局部作用与有限域可限制有效局部对比。没有本轮尺度收敛证据。


23. 本轮不能决定必须改变G/K、界面或损伤模式。单一sc场若无目标，不能据此否定真实组织异质性，也不能同时改多个机制追求成功。


24. 不应立即把不规则删键网络当修复。原始完整非局部图是本轮受控条件；未来如研究不规则离散，须先证明参考m、能量和离散收敛，且不能人工设计目标碎块。


25. 可支持：局部sc确实进入native；均匀回退/历史/能量检查通过；同一开环驱动下的实际功、真实断键、图分量、空间差异及视频。不能超出这一对真实工况。


26. 未标定：300 µm、seed、alpha、q与实际组织组成对应、sc分布及其宏观强度/韧性；当前剪切不是确定的微泡声流；细屑水动力和完整湿润裂面仍未解析。


27. 动态结果图与三个MP4使用真实已接受FP64状态，位移×1。Fig1/2来自冻结初始材料场，Fig3含明确标记的离线规定路径和原生验证，不能把这些输入/离线图称全时程native动力学。每视频152个真实时序帧+60帧末态定格，24 FPS、1920×1080 H.264。各面板显示实际时间，最大两组配对时差11.75 µs；无轨迹插值。


28. U：CPU 2165.470s、wall 2166.185s、峰值RSS 334.406MiB；H：CPU 2138.554s、wall 2139.171s、峰值RSS 336.707MiB。构建、离线路径、静态图及视频成本另见RUN_RESOURCES和各资源JSON；并行wall不相加作为总历时。资源计量覆盖所列计算/后处理进程；依赖安装、交互检查、预览修图及网络传输未全量计时，不能解释为完整服务账单。


## 表面与视频限制

U：末态75/75个自由分量没有可表示的受支持三维湿润表面；未解析节点80；H：末态68/68个自由分量没有可表示的受支持三维湿润表面；未解析节点74。这些分量的流体受力未解析，不能解释为物理上无阻力。当前surface按支持四面体、张开和外部连通及有限遮挡探针近似；未知连续湿润面积无法由点数据确定。该限制两组相同，但会限制后续输运及二次碎裂解释。

视频1比较损伤与真正自由分量；视频2展示近邻存活/软化键与透明失效轨迹，右侧q随原始ID运动；视频3展示真实分量和内部全部存活材料键。颜色由每组精确成员集合注册ID稳定生成，跨U/H同色不表示同一块。近邻显示阈值1.01dx不改变使用全部原始存活键的拓扑；无虚构连续壳、放大位移或视觉聚类。

## 下一轮唯一优先方向

优先验证材料相关尺度在当前非局部离散中的可分辨性：围绕lc/δ≈2及窄试样边界占比，开展受控尺度研究；保持统一起裂与单独sc异质性，先建立足够空间采样与有效局部参数对比证据。暂不同时引入弹性随机场、界面变化或不规则删键网络。

## 文献支持与适用边界


- [Fan et al. 2022 — A Meshfree Peridynamic Model for Brittle Fracture in Randomly Heterogeneous Materials](https://arxiv.org/html/2202.06578v1)：Spatial stochastic material properties and symmetric two-point formulations in heterogeneous LPS. The paper uses harmonic means for elastic parameters and a separate fracture-energy construction. Its brittle critical-stretch law is not the present endpoint-star variational softening law. We borrow spatial-field and symmetric-property principles, not its fracture criterion or convergence claims.


- [Staessens et al. 2020 — Structural analysis of ischemic stroke thrombi: histological indications for therapy resistance](https://haematologica.org/article/view/9534)：Histology of retrieved stroke thrombi identifies interspersed RBC-rich and platelet/fibrin-rich structures with distinct local organization. Histological composition is not mapped to q, sc, or fracture toughness here. No measured patient material field is used.


- [Non contrast enhanced volumetric histology, 2022 — Non contrast enhanced volumetric histology of blood clots through high resolution propagation-based X-ray microtomography](https://www.nature.com/articles/s41598-022-06623-8)：Three-dimensional volumetric clot imaging motivates examining internal spatial organization instead of only bulk composition. Imaging resolution and visible organization do not calibrate the 300 um Gaussian correlation length.


- [Liu et al., JNIS 2022 — Quantification of clot spatial heterogeneity and its impact on thrombectomy](https://pmc.ncbi.nlm.nih.gov/articles/PMC11178127/)：Spatial heterogeneity was quantified using analysis grids spanning 0.2–0.6 mm. Analysis grid size is not a measured Gaussian random-field correlation length. We do not infer that 0.3 mm is a physical correlation length.


- [Eyisoylu et al. 2026 — Arterial shear rate determines the structure and mechanical properties of whole-blood clots](https://pubmed.ncbi.nlm.nih.gov/42576979/)：Local effective stiffness varies with composition; platelet/fibrin-rich regions were stiffer than RBC-rich regions in that study. This supports the existence of local mechanical variation but does not provide our softening thresholds. G and K are deliberately unchanged.


- [Liu et al. 2024 — Adhesive and cohesive fracture of blood clots: Experiments and modeling](https://static1.squarespace.com/static/562c1058e4b0f8ea949c2a94/t/66f453f01994d4758b64b616/1727288309461/1-s2.0-S0022509624003247-main.pdf)：Distinguishes bulk cohesive and interface adhesive fracture. Reports rate-dependent fracture energies and differing RBC/platelet effects on bulk and interface fracture. Neither stiffness nor sc is calibrated Gc. The present kappa-based law has no explicit rate variable; inertia and drag do not reproduce measured intrinsic rate dependence.


- [GSTools 1.7.0 documentation — Gaussian covariance model](https://geostat-framework.readthedocs.io/projects/gstools/en/stable/api/gstools.covmodel.Gaussian.html)：GSTools Gaussian covariance has a rescale convention. This implementation explicitly sets rescale=1, so C_z(r)=exp(-(r/300um)^2). The bounded tanh transform and finite sample alter empirical q correlation. The nominal covariance applies to z.


历史10%异质性实验采用250 µm场、变动s0及球形压入，未产生目标碎块。它说明旧设置的结果，不能作为本轮sc-only剪切的直接均匀对照。

## 交付与复现

`SOURCE_FREEZE.json`、`SOURCE_CHANGES.diff`、`BUILD_AUDIT.json`记录代码/库；`HETEROGENEITY_FIELD.npz`、两个属性CSV和场元数据记录完整实现；两个NATIVE_RESULT、ENERGY_LEDGER、NATIVE_EVENT_HISTORY与TRUE_FRAGMENT_HISTORY记录接受证据；`REPRODUCE.md`给出环境和执行命令；`OPEN_RESULTS.html`链接图、视频与源数据。历史数据未删除；没有优化、微泡CFD或完整V3生产运行。


绘图曾因空失效键集合报错，VTK视频曾在重建坐标控件后崩溃；修复仅限后处理并保留失败日志，没有重跑原生工况。失败渲染及短帧渲染检查未计入成功后处理资源总计。
