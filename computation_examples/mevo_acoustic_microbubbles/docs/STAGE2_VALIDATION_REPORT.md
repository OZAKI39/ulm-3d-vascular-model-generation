# Stage2 validation report — 2026-10-01

## 1. Executive summary

本轮已在现有工程和服务器中实际实现、编译、运行并测试。**软件链路有显著进展，Stage2整体验收尚未通过。**
APECSS非线性脂质壳后端、响应表、单泡平移参考、N=2双泡参考、适用范围保护与CPU/GPU比较均已落地。
Zhang2011的吸引/排斥及共振峰趋势得到部分复现，但图3的间距导致力符号反转未复现。
按用户指定闸门，Lazarus三档聚集、圆管验证演示与五频扫频**没有运行**；不将旧线性替代模型结果改称新物理结果。

最终全套测试34 passed；响应表/validity修改后的相关15项与黏度映射专项检查另存日志。
ruff检查通过。测试通过仅表示这些软件/参考检查通过，不等于实验校准或近壁物理准确。

| 阶段 | 实际状态 |
|---|---|
|2A|审计、基线tag、逐文件核对完成|
|2B|APECSS官方构建/算例/参考数组通过|
|2C|PROTEUS源码/许可证/复用审计完成；未声称执行MATLAB|
|2D|SonoVue参数出处与有界近似数字化完成|
|2E|非线性径向软件验证通过；表格仅部分单元可用|
|2F|单泡参考与平均力数值检查通过；准稳态误差明确|
|2G|部分验证；Zhang图3失败保留|
|2H|全文提取完成；50/75/100kPa模拟因闸门未运行|
|2I|实际CPU/GPU性能和交叉验证完成；C++依赖/适配修复|
|2J|因核心验证未通过而未运行；没有新增频率基底|

## 2. Existing code retained

基线提交`fc3fac0e6fe7f0a57bd90812800714a97fff94a5`保留为`stage1-demo-baseline`，开发分支`stage2-validated-physics`。
保留PZT控制、基底HDF5/复相位约定、有界声场插值、FlowField.velocity(x,t)、已有几何接触、VTK输出与原始结果。
声学修改仅针对已确认的0.6.2 C++接口问题；CuPy求解过程和缓存身份未改变。

实际回归：**18,377个原始结果文件零缺失、零哈希变化**；接手时全交付清单18,440文件也已核对。
缓存基底重新组合对比已保存direct场：1MHz误差9.78255119e-8，2MHz误差1.56248854e-7，与基线一致。
150泡迁移10ms重跑：所有已保存时刻一致，全部位置最大坐标差**0m**。
两个声场网格重新执行共同坐标比较，原有压力/梯度差异保留；没有重复计算基底。
证据：[REGRESSION_STATE.json](../results/stage2/audit/REGRESSION_STATE.json)。

## 3. Third-party software audited

APECSS v1.8/commit`a0ee071ca247b5489e997d90e83b5f53c40e4ec7`；
PROTEUS commit`fc24fcbcf9623006e8f79989c0126c45a76bc23b`；
实际k-wave-python0.6.2、CuPy14.2.0、NumPy2.2.6、SciPy1.15.3。
官方实现/作者代码/补充材料先于Level3参考实现调查。
Doinikov2015作者Mathematica补充代码确实存在，但下载不可用；没有把这一情况写成“不存在代码”。
完整记录：[THIRD_PARTY.md](../third_party/THIRD_PARTY.md)、[许可证](../THIRD_PARTY_LICENSES.md)、[REFERENCES.bib](../REFERENCES.bib)。

## 4. APECSS validation

实际gcc13.3.0，Release、O3、march=native、double8字节，无fast-math。
7个ultrasound算例加Rayleigh算例及8个官方绘图脚本完成。
4份官方参考数组均满足上游rtol1e-6/atol3e-16，relativeL2为2.74e-14至4.05e-12。
未修改上游源码。正式通过记录：[APECSS STATE](../results/stage2/apecss_validation/STATE.json)。
早期源文件同步/示例文件名错误保留，不纳入通过数量。

脂质模型是APECSS可逆三段Marmottant特例，破裂阈值取水表面张力；不是2005论文所有破裂历史分支。
RPAR气体声辐射项和壳黏性逐项审计；不另写通用Python径向库。
上游驱动为−Pa sinωt，对应+iPa相量；已正确转换到项目cos相位并测试。
[模型审计](APECSS_MODEL_AUDIT.md)。

## 5. PROTEUS reuse decision

保留其MIT和内嵌BSD/MPL通知。非线性响应、Segers参数、flowode23、稀疏网格/VTU处理及完整CEUS链均已阅读。
其MATLAB流程并不提供本项目所需的通用四面体Python插值或本轮两个耦合参考方程。
继续复用APECSS和已有PyVista接口，未整包并入PROTEUS，未复制其求解源码。
[组件决策表](PROTEUS_REUSE_AUDIT.md)。

## 6. SonoVue parameter provenance

Tu2009 Table2：R1.78µm、χ0.30N/m、κ3.2e-9kg/s、150kPa/2.5MHz。
Table3是同一观测的替代拟合：R1.72µm、χ0.22、κ2.5e-9。
**40kPa/1.70µm/Table1是BR14，不能改称SonoVue。**
Fig4选点/像素标定/图片SHA和线性χ、logκ插值可追溯；拒绝超数字化半径范围外插。
±0.02N/m和±0.06logκ是数字化分辨率假设，不是生物群体标准差。
σinit=.036N/m是假设，另运行0/.072敏感性算例。
Tu的有效径向黏度.002Pa·s包含热阻尼；新增显式输入可与平移阻力所用物理黏度分开。
[参数文档](SONOVUE_PARAMETERS.md)。

## 7. Single-bubble validation

正式粗表96个任务全部完成，60节点中59个周期/弱压缩有效，24个插值单元**无一通过**5%中点检验。
聚焦表201个任务全部完成：R1.6–1.96µm，f2.3–2.7MHz，p5–40kPa；125节点有效，**39/64单元可查询**。
其余单元明确拒绝；没有把完整数值范围写成完整物理适用范围。
选定细步/长时间对照的最大谐波L2变化约1.41e-7；包装器与官方lipid波形L2为1.34e-8。
不存在实验SonoVue曲线拟合或完整诊断脉冲复现的声明。

[Figure A](../results/stage2/single_bubble/20261001T201550430625Z/Figure_A_radial.png)；
[响应表状态](../results/stage2/single_bubble/STATE.json)。
直接外部ODE仅离线执行；慢时间轨迹批量查询复谐波表，无逐泡MHz积分。

## 8. Translation validation

Doinikov–Dayton2006厚黏弹壳N=1参考完成22个粒径/压力/频率/壳厚/细步任务。
弱场平均力对Eq63误差为**0.343%、0.196%、0.057%**。
20周期、R2µm/150kPa：完整模型位移0.0970µm，使用相同周期力的准稳态模型0.1025µm，约5.64%差异；
R3µm/150kPa约12.4%差异，且最大Re≈1.71。低压力暂态与局部正负位移保持原值。
细步最大位移差3.27e-21m表示当前自适应容差下数值稳定，不能代表模型准确到该尺度。

单频梯度下，非线性V的基频投影与完整周期平均力是数学等价关系，实测最大相对差5.27e-13。
多谐波梯度检查另验证谐波求和；没有人为制造“非线性导致单频公式失效”。
[Figure B](../results/stage2/translation/20261001T202344093431Z/Figure_B_translation.png)、[状态](../results/stage2/translation/STATE.json)。
2026交叉核对确认为Physical Review Fluids11,033606；其history-force/thermal模型不被自动替换进当前生产模型。

## 9. Two-bubble validation

Zhang：完成21个对照任务。R1=4µm/R2=6µm吸引，3µm/6µm排斥；半径扫描峰位2.9µm，与代表性趋势一致。
64周期步长对照最大间距差4.30e-16m。然而R3/6µm在120µm时fB仍为−6.14e-18N·m²，
论文图3显示正值，**此项失败**。只计算图4前0.2ms，未复现最终稳定泡对。

Doinikov2015：4个弱场与4个近场阶数任务；2.15/2.95MHz运动方向符合图3，M8/M12接触前最大间距差2.37e-10m。
Gaussian宽度假设、缩短弱场脉冲、无包膜/形变/润滑等限制均记录。接触终止而不宣称接触后预测有效。
总体为**partially validated**，不能据此解除H/J闸门。

[Figure C 双泡参考](../results/stage2/pair_zhang2011/20261001T203554552347Z/Figure_C_pair.png)、
[独立响应远场对照](../results/stage2/pair_zhang2011/20261001T203554552347Z/Figure_C_reference_vs_farfield.png)、
[Zhang状态](../results/stage2/pair_zhang2011/STATE.json)、[Doinikov状态](../results/stage2/pair_doinikov2015/STATE.json)。

## 10. Clustering validation

Lazarus全文含图像公式已读。来源是近壁二维自制PFB微泡，固定平均半径/压缩率β=5e-7Pa^-1，
时均吸引力加Stokes阻力和附加质量；100nm接触、0.5ms经验重排，不是非线性SonoVue壳求解。
其论文明确在远场公式适用距离外继续使用该式。已提取接触、体积加权位置、A(1−exp(−Bt))和sqrt(S/N)归一化定义。
50/75/100kPa均保存**NOT_RUN_STAGE2G_GATE**。除闸门外，还缺实测群体坐标/半径及可校准壁面/历史阻力。

[全文提取与原因](../validation/lazarus2017_clustering/README.md)、
[Figure D 未运行状态](../results/stage2/lazarus2017/Figure_D_gate_status.png)。
没有新聚集率、聚集尺寸或压力曲线；旧32泡几何cluster不能填入此处。

## 11. Numerical convergence

已做径向采样/稳态时长、单泡平移步长、Zhang步长与Doinikov多极阶数检查。
Stage1共同坐标压力网格差仍为迁移3.43%、聚集1.40%；原聚集20ms轨迹dt比较末位置RMS差5.72µm仍成立。
局部参数表验证和全局守恒/CPU-GPU一致性不替代局部物理或实验验证。
全部测试、初期失败及后续修正日志保存在`results/stage2/audit/`。

## 12. Performance

以下为8谐波、无接触合成阵列的中位耗时，不是聚集轨迹。GPU整步含CPU接触及传输。

|N|NumPy力ms|CuPy力ms|NumPy合成步ms|CuPy混合步ms|
|---:|---:|---:|---:|---:|
|32|0.044|0.435|0.127|0.582|
|64|0.084|0.417|0.190|0.619|
|128|0.213|0.425|0.446|0.725|
|256|0.838|0.429|1.531|1.240|
|512|2.760|0.459|6.321|4.117|
|1000|10.026|0.728|19.439|9.417|

CPU/GPU最大力相对差3.58e-16。N1000显存池保留量约108MB，不是精确VRAM峰值。
密集接触链N32/64/128的CPU约0.626/2.112/9.824ms，证明接触可主导成本。
万泡响应查询NumPy4.681ms/CuPy3.871ms，同点误差0；这是有限吞吐测量，不是全表误差证明。
保留小N NumPy参考；当前N≤1000可用有上限GPU全对。默认上限2000；不实现未经物理误差评估的力截断。
扩大规模前优先评估接触候选邻域搜索；本轮不凭猜测编写复杂邻居算法。
[性能原始表](../results/stage2/performance/STATE.json)。

## 13. GPU/CPU utilization

服务器Ryzen7800X3D，8物理/16逻辑核，1NUMA，cgroup约7.68核，无Slurm；APECSS池选择6worker、每worker1线程。
聚焦表11.45s/201任务，17.56case/s，平均3.10核，采样聚合RSS约1.02GB。约52%worker利用率，未宣称CPU满载。
RTX4090实际执行声场、pair力和响应插值。C++两档小网格对比误差3.26e-7/4.37e-7，
单次含I/O耗时约为CuPy的56%–58%；生产大网格尚未比较，CuPy保持默认。
[后端修复/benchmark](KWAVE_BACKEND_BENCHMARK.md)、[硬件记录](../results/stage2/system/SERVER_DOCTOR.json)。

## 14. Known limitations

核心限制包括Zhang图3失败、响应表空洞、壳初始张力/数字化不确定性、有效径向黏度映射、
无history drag、有限Re、无近壁动力学、N=2不可直接推广N体，以及原密集接触轨迹未收敛。
完整列表与保留错误：[KNOWN_ISSUES.md](../KNOWN_ISSUES.md)。无隐藏epsilon修补、力放大或为图像调参。

## 15. What is now safe to use

可以使用现有声学缓存进行兼容控制组合；运行APECSS官方基准；进行有出处的参数敏感性研究；
在声明假设及经过检查的表格区域内进行软件/方法开发；使用strict拒绝越界、保存状态；
使用受限远场CPU/GPU核、文献参考对照和VTU空间/时间插值测试。
这里的“可用”限定于明确模型假设和软件验证范围，不是临床或患者特异预测资格。

## 16. What is still experimental

SonoVue定量预测、近场耦合/多泡密集聚集、壁面作用、Lazarus重现、圆管声控效果和五频扫频仍为实验性或未运行。
LegacyRegularizedBjerknes及所有旧linear_surrogate配置明确demo only；原Marmottant占位类继续显式拒绝调用。
没有不必要地重写完整MarmottantReferenceModel，因为APECSS已提供本轮选择的成熟可逆特例。
验收A–J、M–R中相应软件/出处/接口项已有证据；K仅代表性趋势部分通过，L记录清楚未复现原因。
因此**不宣称本轮所有科学验收条件已满足**。

## 17. Readiness for SimVascular

FlowField.velocity(x,t)保持；ZeroFlow、Poiseuille、原静态SimVascularFlowField仍通过。
新增通用VTKFlowField支持可配置点速度、坐标/速度SI缩放、时序线性插值与域外拒绝，synthetic VTU测试通过。
WallInteractionModel、GeometricNoPenetration和显式TODO的DoinikovWallReference接口已建立。
当前达到**接口准备完成、患者物理预测未准备完成**的状态。没有启动患者MeVO或新的BraVa轨迹，
没有恢复旧小鼠WSS/压力敏感性/微泡批次，也没有修改迁移兼容链接。
下一步首先独立核对Zhang公式实现与图3所用参数/平均窗口，并进一步细化响应表；通过闸门后再进入H/J。

复现入口：[STAGE2_QUICKSTART.md](STAGE2_QUICKSTART.md)。源模型追溯：[MODEL_MAP.md](MODEL_MAP.md)、[PHYSICS_PROVENANCE.md](PHYSICS_PROVENANCE.md)。
