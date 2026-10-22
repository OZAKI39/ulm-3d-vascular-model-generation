# s1-11 实际承托闭合与安全拆除检查

状态：**NEEDS_ADJUSTMENT**。生产候选：无。 当前文件仅用于诊断。

1. 中心线实际承托是否 100%：否，98.74371859%。
2. 延长端口是否 100%：否，99.17355372%。
3. 完整逐层区域是否 100%：否，99.36224476%。
4. 最长未通过承托检查的弧长：0.21980531 mm；失败采样点 25 个。过近顶面也作为失败，不能用较低幕墙层代替。
5. 支撑自身未通过上一层承托检查的区域：284 个。

## 局部修补与停止条件

本轮读取既有报告、逐点覆盖、多边形、拆除及空移记录；保留原姿态搜索结果，只处理 O002 与既有两个候选姿态。完成 32 次有真实切片指标的修补。每条路线最多 8 次；仅在连续两次缺点数和缺口面积均不改善时提前停止。运行异常单独记录，不算作收敛。

所有变化仅发生在支撑的实际缺口附近，尺寸采用当前线宽或层高的倍数，网格顶部的半层偏置来自切片平面位置；没有放宽覆盖阈值或改动血管。缺口数据和每次变化见 `support_closure_iterations.csv` 及各候选的 `local_repairs.json`。

本次诊断展示选择 `O000_R20_R3_C3_LAYER_I05`，排序首先比较三类覆盖率中最低的一项。它只用于呈现失败证据，不是已验收的生产候选；局部迭代是否接受仍要求缺点数与缺口面积均不退步。

覆盖未达到三个 100% 且支撑自身未闭合时，不开始修改拆除几何。下面仅为当前诊断模型的只读虚拟剪切、载荷路径与顺序检查；没有将其冒充生产拆卸方案。

6. 计划支撑段数：35。
7. 有限虚拟剪切后的实际组件数：1。
8. 最大组件接触弧长：216.515 mm；包围盒各轴上限：[97.69075012207031, 72.99220275878906, 20.769990921020508] mm。若仍为巨大连通组件，则不满足安全拆卸要求。
9. 候选剪切 20 次，其中几何安全检查通过 0 次；不是已批准操作次数。
10. 不要求所有剪口初始可达。每次剪切前保留尚未切除的材料，每次移除前保留其它组件，随后更新障碍物。

## 顺序与血管载荷风险

- 第 1 轮：新增可达剪口 []；移除 []；剩余 ['G000']。

11. 不允许通过血管弯曲来获得移除路径。未找到刚体路径表示“未能证明可拆”，不等于已证明必须弯曲血管。
12. 检出的闭合截面套环：0 个。
13. 最大局部包覆角：225.687°。
14. 超过 120° 的组件：1 个。
15. 超过 150° 的组件：1 个，默认失败；本轮未启用人工审核例外。
16. 分叉附近当前接触组件最大弧长：216.515 mm。
17. 最大几何力臂：57.200 mm。
18. Weak-link ratio：尚未优化或通过真实切片确认，保留为空。实际血管截面由冻结 STL 截取；不使用半径推算，也不将面积比解释为断裂力比。

## 验证边界

源保护核对 485 个文件，改变 0 个。真实 Bambu 切片成功：True；自动支撑关闭：True。未解决关键警告：['It seems object BG001 vascular core + custom support has floating regions. Please re-orient the object or enable support generation.']。

软件只能检查几何、实际挤出和碰撞代理，不能保证 ABS 一定先断或血管一定不会断。当前不建议完整实物打印。拆除原则仍为先剪断、再移除长幕墙、最后处理短接触块；未认证的位置不能直接横向掰动。

## 未通过的硬门槛

- COVERAGE_100
- SUPPORT_SELF_ZERO
- NO_VESSEL_DEFLECTION_REQUIRED
- COMPLETE_SEQUENTIAL_REMOVAL
- NO_LONG_LEVER_CONTACT
- WRAP_AT_MOST_150
- HIGH_WRAP_EXCEPTION_JUSTIFIED
- WEAK_LINKS_VALIDATED
- WARNINGS_RESOLVED

## 文件入口

`selected/`：血管 STL、自定义支撑 STL、仅拼接的审阅 STL 和原生 `X1C_project.3mf`。血管与支撑保持两个 normal part。

`QC_SAFE/`：17 张图。图 03 文件名沿用任务要求，但只要失败点存在就显示红点及真实百分比，不能画成全绿。

`safe_removal_sequence.json`：生产方案为空；诊断检查与未完成原因单独记录。

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-11_vascular_core_support.py --stage closure
```

## 基线对照与真实运行证据

原基线：`O002_R25_C2_LAYER_BASE`。各项独立列出，不能以逐层覆盖改善替代中心线及支撑自身合格。

|指标|本轮开始时的基线|最终诊断候选|
|---|---:|---:|
|中心线覆盖 / %|98.89447236|98.74371859|
|端口覆盖 / %|99.48347107|99.17355372|
|逐层覆盖 / %|98.42081545|99.36224476|
|失败采样点|22|25|
|逐层缺口面积累加 / mm²|1.32643219|0.50847367|

停止记录（运行异常与真实停滞分开）：

|路线|轮次|原因|
|---|---:|---|
|O002_R25_C2_LAYER_BASE|2|TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT|
|O002_R20_C2|7|TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT|
|O000_R20_R3_C3_LAYER|1|LOCAL_REPAIR_RUNTIME_FAILURE_NEXT_LOCAL_VARIANT：SUPPORT_ASCII_PRECISION_ROUNDTRIP_FAILED|
|O000_R20_R3_C3_LAYER|2|LOCAL_REPAIR_RUNTIME_FAILURE_NEXT_LOCAL_VARIANT：SUPPORT_ASCII_PRECISION_ROUNDTRIP_FAILED|
|O000_R20_R3_C3_LAYER|3|LOCAL_REPAIR_RUNTIME_FAILURE_NEXT_LOCAL_VARIANT：SUPPORT_ASCII_PRECISION_ROUNDTRIP_FAILED|
|O000_R20_R3_C3_LAYER|7|TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT|
|O000_R25_C2|3|LOCAL_REPAIR_RUNTIME_FAILURE_NEXT_LOCAL_VARIANT：SUPPORT_ASCII_PRECISION_ROUNDTRIP_FAILED|
|O000_R25_C2|4|LOCAL_REPAIR_RUNTIME_FAILURE_NEXT_LOCAL_VARIANT：SUPPORT_ASCII_PRECISION_ROUNDTRIP_FAILED|
|O000_R25_C2|7|TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT|
|O001_R20_C2|7|TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT|
|O001_R25_C2|7|TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT|

局部修补实际接受 0 次，未通过回退条件的试验全部保留供追踪。当前诊断的 25 个失败点中，22 个存在过近的顶面，其余 3 个缺少有效顶面承托；不能把这两类都描述成完全裸悬空。

支撑自身悬空区域：原基线 45 个，当前诊断 284 个。任何单项改善均不能抵消这一未通过状态。

永久测试：199 项，失败 0，错误 0，跳过 0。测试通过只说明检查逻辑按预期执行，不代表制造门槛通过。

真实合成连接片 `actual_thin_web`：连接片保留判据 FAIL。详细逐层证据见 `validation/actual_thin_web/`；连接片保留不能替代完整覆盖及支撑自身检查。

真实合成连接片 `actual_thin_web_C3`：连接片保留判据 PASS。详细逐层证据见 `validation/actual_thin_web_C3/`；连接片保留不能替代完整覆盖及支撑自身检查。

Bambu 的对象级悬浮警告通过候选目录中的 `warning_localization.json` 对应到实际失败采样点、逐层缺口及支撑自身悬空区域；这些区域未闭合前均保留 UNRESOLVED。

实际刀路空移审计：78251 段空移，12050 段触发近碰代理。结果使用逐时刻已打印几何及喷头代理，属于几何风险筛查，不是实测碰撞或已校准的失效概率。明细见所选候选的 `travel_near_miss.json`。
