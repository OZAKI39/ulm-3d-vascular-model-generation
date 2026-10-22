# BCF V1 历史说明：独立血管芯与无底浇注夹具

> 本页描述已退休的穿墙端口/U 槽版本，保留用于理解旧输出。当前 `s1-9` 命令默认运行 V2，请使用 [BCF V2 内墙支撑与切割说明](BCF_V2内墙支撑与切割使用说明.md)，不要将下文旧打印文件作为当前生产件。

本轮生产路线改为：**先单独打印并清理 ABS 血管芯，再从顶部装入无底 PLA 框架；培养皿光滑底面直接形成 PDMS 底面。** 旧 V1/V2/V3/V3.1/RFUS 支撑研究不再用于生产选择。现有全脑、ROI 可视化入口保持原样。

## 使用入口

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-9_vascular_casting_fixture.py
```

此命令读取冻结生产源、生成夹具 CAD、检查装配、公差与接触、在本机 Bambu 实际切片，再生成报告和 20 张 QC 图。它不重跑语义迁移、ROI 选择、半径补偿、VascularMD 或端口光顺，不控制打印机。环境为本项目 `run_python.sh` 指向的共享 Python；依赖已在 requirements-printing.txt、requirements-fixture.txt 中固定。

只检查几何可加 `--geometry-only`，此时没有本轮切片证据，状态保持 NEEDS_ADJUSTMENT。`--config` 指定 YAML；`--output` 指定另一输出目录。`--no-step` 和 `--no-render` 只用于诊断，正式交付不要省略这些产物。

## 输入与保护

自动发现 s1-7 当前已选 `candidate_0008` 的 `manufacturing/core_with_ports.stl`，逐一核对上游文件哈希。原血管在 `vascular_core/BG001_RMCA_core_with_extensions.stl` 中为字节相同副本；位置变换单独写入 `assembly/casting_core_transform.json`，缩放系数始终为 1。

框架内腔根据冻结中央血管表面、完整端口范围和直段可夹持区重新确定。当前约 81.765 × 69.976 × 41.011 mm，侧壁厚 2 mm、向外法兰宽 5 mm、底面血管间距 10 mm、顶部净距 12 mm。顶部及底部都开放。

端口数量从冻结 `geometry.json` 读取。保守使用 target → end 的已知直段，在夹持区作三个真实 STL 截面复核。当前 O2 为主基准，外侧端面挡块提供轴向参考；I1、O1、O3 的贯通槽保留 ±1 mm 轴向余量。径向间隙 0.22 mm 来自实际 0.22 mm 线宽的候选。59 组确定性误差下，主基准加浮动方案通过 59 组，全固定方案通过 17 组；误差包络是工程假设，尚未经实物校准。

## 文件去向

以下路径都相对于：

`outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_casting_fixture/`

|目录或文件|内容与用途|
|---|---|
|vascular_core/|独立 ABS 血管芯原样副本|
|frame/bottomless_PLA_casting_frame.stl、.step|唯一 PLA 无底框架|
|closures/port_*_closure.stl、.step|4 个可拆闭合插片；O2 插片带外侧基准挡块|
|closures/port_*_locking_pin.stl、.step|4 个外侧锁销，限制插片被向上抬起|
|closures/closures_and_pins_print_layout.stl|按已验证朝向摆放的插片与锁销，零件之间保留间距|
|assembly/full_fixture_assembled.step|独立零件装配；血管保留原三角面的 faceted BRep，不重新拟合|
|assembly/full_fixture_review.stl|仅供审阅的多零件 STL，不应当作同一种材料的一体打印件|
|assembly/assembly_sequence.json|从顶部放芯、先主基准再其它插片及锁销的顺序|
|bambu/frame/sliced.3mf|真实 PLA 框架切片，Auto Support OFF|
|bambu/closures/sliced.3mf|真实 PLA 插片和锁销切片，Auto Support OFF|
|bambu/ABS_core_reference/sliced.3mf|独立 ABS 芯参考切片，沿用现有自动支撑工艺|
|config/clearance_study/|9 个公差组合、59 组确定性误差、另两个径向间隙的 CAD 候选|
|tables/|端口、装配、公差、接触、间距、密封及插入检查表|
|QC/|20 张模型、无底剖视、装配、接触和爆炸图|
|reports/|永久测试、最终保护审计及本轮命令日志|
|vascular_casting_fixture_summary.json|机器可读结论与失败门槛|
|vascular_casting_fixture_report.md|中文结果和实际制造顺序|

Bambu 的几何输入 STL 与打印坐标变换也留在各切片目录。无喷嘴/材料切换的隐式操作；框架和插片使用本机 X1 Carbon、0.2 mm 喷嘴的 PLA 模板，血管参考片使用 ABS Vessel。没有上传、发送或开始打印。

## 培养皿与实物边界

本轮没有找到实际培养皿尺寸。完整法兰需要至少 **127.370 mm 的平坦圆形底面直径**；这不是培养皿标称外径。测量后可新增 `config/casting_dish_measurements.json`，使用数值键 `dish_inner_flat_diameter`、`dish_wall_height`、`dish_floor_thickness`（mm），再次运行进行尺寸核验。不要填猜测值。

下托座、插片和锁销仅在侧壁端口直段附近限位。几何检查中，中央血管没有夹具接触，内腔无非法 PLA 侵入。名义装配有间隙，4 个捕获位置并非 4 处预压接触；重力或浮力作用后可在有限间隙内到达限位面。

外侧端口环槽、U 槽周边槽和法兰底部闭合槽是密封几何，没有认证胶黏剂。真实密封、ABS/PLA 收缩、插片与锁销强度、抗浮力和拆卸过程仍待实物试验。培养皿未知时的软件状态是 READY_FOR_HUMAN_REVIEW_WITH_DISH_FIT_PENDING，不能写为制造认证或防漏认证。

## 代码职责与测试

入口 `s1-9_vascular_casting_fixture.py` 调用 `vascular_processing/casting_fixture_pipeline.py`。geometry 模块负责冻结输入摆放与 CadQuery 参数化零件；audit 模块负责 FCL 整路径检查、接触、内腔侵入与刚体限位；profiles 模块只做实际本机预设读取与本地切片；tolerance 模块比较确定性公差；review 模块生成独立 QC 和报告；synthetic 模块生成十类永久回归数据。

```bash
./run_python.sh -m pytest -q tests/test_s1_9_vascular_casting_fixture.py
```

清理账本在根目录 `cleanup_inventory.json`、`cleanup_plan.md`、`cleanup_result.json`；阶段回顾在 `research_retrospective/SUPPORT_RND_RETROSPECTIVE.md`。本项目没有 Git 工作树，旧算法源码因此保留，旧支撑候选几何、切片和缓存已删除；历史报告中这些大文件的旧链接不再是现存交付物。
