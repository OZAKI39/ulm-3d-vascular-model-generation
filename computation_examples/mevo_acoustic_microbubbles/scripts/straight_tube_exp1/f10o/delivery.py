"""Write concise reports from archived quantities, without extending scientific claims."""

from .common import BASE, ROOT, read, write, sha


def main():
    g = read(BASE / "analysis/FINAL_GATE.json")
    c = read(BASE / "analysis/COMPARISONS.json")
    runs = read(BASE / "analysis/RUN_ANALYSIS.json")
    w = runs["validation_v1/WEAK_DT1"]
    matched = read(BASE / "analysis/MATCHED_0p1MS.json")["data"]
    commits = {
        name: read(BASE / name / "CAMPAIGN.json")["identity"]["execution_commit"]
        for name in ["pilot_v1", "validation_v1"]
    }
    text = f"""# F10O：清理、CPU B2复核与有条件XY控制

已完成安全清理、既有B2接入、有限科研计算、失败审查和四组英文PNG/PDF。**强控制未通过原重叠预算，未启动XY优化或长窗口。**这是执行了预登记停止规则的完整交付，不是导航成功。保留EXPLORATORY、PHYSICS_UNCALIBRATED、training_ready=false。

## 基线与执行身份

完整F10N基线e2dfb4f41cea33082583bbf8424a963ee194e85c，独立分支dev/fast-pzt-f10o-b2-xy-control-20261013；原F10N工作区不变。

- 独立清理提交：19c11e885。
- B2工程接入提交：3f5a117c1。
- 150步试运行提交：{commits["pilot_v1"]}。
- 正式有限数值对照提交：{commits["validation_v1"]}。
- 后处理、绘图及最终文档有自己的Git身份，不作为科学积分版本。

服务器/workspace/mevo_acoustic_microbubbles/f10o_b2_control，Supervisor独立标签mevo_f10o_pilot_v1和mevo_f10o_validation_v1，均有330 s外层超时及150 s单条预算。源代码Git归档不可变，全部输入在运行前核对SHA；四路控制直接读取F10N INPUT.json完整精度，未根据报告舍入值重建。

## 清理结果与保护

实际删除src/mevo_mb/fast_sim/vessel_kernels.py中的无调用f64转换函数，共5行。AST名字计数、模块消费者、配置、入口、动态导入/属性搜索以及删除后原F10N完整只读验证共同支持删除；服务器还实际编译了清理后的Warp模块。没有把无人导入的CLI当作无效文件。清理幅度很小，主要通过docs/F10O_ACTIVE_CODE_MAP.md减少开发上下文，不宣称大量历史代码被淘汰。

CPU B2、GPU诊断、R3真实壁面、旧Y管/三PZT继承依赖、历史命令入口和拒绝未实现物理的测试桩全部保留；报告、HDF5、场、响应表、失败和哈希均保留。详见F10O_DEPRECATED_CODE_AUDIT.md及cleanup/原始diff。没有改动任何历史执行工作区。

## 实际模型与复用

32固定ID、R=1.78 µm、2.5 MHz、四枚7 mm源、均匀水圆片复压力和一致梯度、原APECSS及F2近场策略不变。缓存3/1.5 µm直接复用，无k-Wave新求解，无CFD。腔体充水、背景流严格零，双顶部孔及两端封闭几何仍在；此为控制代理，不是实验校准。

WaterB2每步调用原F10H geometric_contact和F10G resistance。刚度保持0.04641319518378658 N/m，非线性残差1e-13 m、最多24次更新、20次线搜索，都是F10K轨迹参考原值。没有重写PCG、本构或润滑。GPU仍计算原声学力；同一R3 cavity_project对预测及B2提案均检查。任何非零壁修正直接拒绝这个远壁参考，不串接PBD。14项实际接口核对通过，包括CPU/GPU三种控制力输入、封闭X拒绝、顶部数值域拒绝、低于/高于APECSS表拒绝。

M0采用两次局部润滑迭代及12次PBD投影；B2求解原成对阻力矩阵并采用几何隐式单边弹簧。因此接触和迁移率近似一起改变，跨模型差异不能称为纯数值误差，也不能仅归因于弹簧替代。没有据此证明B2更接近真实实验。

## 计算成本及强控制失败

先各150步：OFF/弱/强平均每步分别约3.340/3.637/5.843 ms；九组2 ms完整对照预计102.55 s，低于首轮300 s预算才启动。实际正式OFF DT1/弱DT1/强DT1平均每步分别为{runs["validation_v1/OFF_DT1"]["mean_step_s"] * 1e3:.3f}/{w["mean_step_s"] * 1e3:.3f}/{runs["validation_v1/STRONG_DT1"]["mean_step_s"] * 1e3:.3f} ms；强控制的正式值只覆盖拒绝前短窗。

强控制dt=1 µs、0.5 µs和细缓存三组，在0.175 ms提案处分别达到17.805391、17.820427、17.805390 nm重叠，超过不可放宽的17.8 nm（0.01R）。数学残差分别约3.06e-15、5.94e-16、3.06e-15 m，通过原停止准则；壁修正为0，单次接触迭代修正也低于0.1R。失败是软接触允许的几何重叠超预算，不是线性求解不收敛或墙面穿透。越限虽小仍拒绝，不能改刚度、容差或接纳提案来制造PASS。失败提案、最后合法状态和全部逐步量均已保存。

有效强控制前缀分别为0.174、0.1745、0.174 ms；dt1前缀质心位移约(0.0195296,0.00355137,0) µm。不能称完成2 ms强控制，更不能把此速度外推到30 ms。

## 同ID、缓存与模型差异

共同合法0–0.174 ms前缀中，强控制最大同ID步长差{c["STRONG_DT05"]["fixed_id_m"] * 1e9:.8f} nm，缓存差{c["STRONG_FINE"]["fixed_id_m"] * 1e9:.8f} nm，远低于50 nm门槛，但**完整2 ms准入仍失败**。质心、Kabsch、Rg、张量和比例差独立保存。

在双方都实际保存的0.1 ms同一时刻，M0同ID步长差{matched["M0"]["dt_same_ID_m"] * 1e9:.6f} nm，B2为{matched["B2"]["dt_same_ID_m"] * 1e9:.6f} nm，说明短前缀步长敏感性减小；这不能与F10N完整2 ms的6.416 µm直接作同窗比较，也没有解决完整强轨迹。缓存差同窗约0.000410/0.000402 nm，没有明显差异，长窗未知。

M0/B2强控制只能比较共同保存至0.1 ms：同ID差15.9549 nm、Rg差13.2558 nm，COM差约1.24e-14 m；这是不同接触/迁移率模型的对照。连通性本轮均保持一个分量；B2所有有效前缀无瞬时边事件。OFF/弱控制网络维持完整2 ms，强控制观察不足2 ms；不能把“无事件”写成新持续聚集或物理融合。

## 弱控制形态复核

OFF三组均完成2 ms并严格静止。弱控制三组均完成2 ms且通过原数学、重叠及壁预算；最大同ID步长差{c["WEAK_DT05"]["fixed_id_m"] * 1e9:.8f} nm，缓存差{c["WEAK_FINE"]["fixed_id_m"] * 1e9:.8f} nm。弱控制COM约(0.00160696,0.000197992,0) µm。

Rg从5.903887306降至{w["rg_initial_final_m"][1] * 1e6:.9f} µm，减少{-g["weak_Rg_change_m"] * 1e9:.6f} nm；原M0增加11.1069 nm。B2 X/Y比由1升至{w["xy_initial_final"][1]:.12f}（约{g["weak_xy_change"] * 100:.7f}%），原M0约0.4245%；全局长短比由2.236067977至{w["aspect_initial_final"][1]:.9f}，终点主轴距X约42.9933°。初始X/Y简并，不能解释成从一个确定初始主轴旋转。

这些B2小变化超过本次步长/缓存差异，低成本形态复核有数值证据；但1.25目标未达到、15°方向要求也未达到。此复核复用第二阶段必要弱控制对照，没有在失败之后额外开展形态优化。APECSS表插值、近场声力、接触与壁面水动力仍未校准，有限步长一致不等于完整物理误差收敛。

## XY目标与停止

运行前已登记X=(0.2,0,0) µm、Y=(0,0.2,0) µm，容差0.075 µm；原F10N (3,2,1) µm完整保留。强控制质量门槛失败，因此新的BASELINE/OPTIMIZED XY策略、分段反馈以及5/10/30 ms均为NOT_RUN。OFF和弱控制的2 ms数据是必要数值对照，强控制合法前缀的朝X移动不是新任务导航成功；也没有证明可独立控制X/Y，更未证明全局Z不可控。

数值对照中的X目标误差确有短暂减少：强控制dt1有效前缀约0.200→0.180505 µm，仍大于容差且轨迹被拒；弱控制2 ms约0.198393 µm。这不能支持持续导航或长时间可达性。

## 成本、图件与归档

两个实际数值进程计时合计{g["numerical_process_total_s"]:.6f} s，其中逐轨迹积分及逐步证据{g["integration_s"]:.6f} s；CPU B2调用合计{g["CPU_solver_sum_s"]:.6f} s、阻力组装{g["CPU_resistance_sum_s"]:.6f} s。GPU用于声场/力和壁检查，未单独进行核计时，不能把主机余量称为GPU纯计算时间。峰值进程RSS {g["peak_process_rss_MiB"]:.3f} MiB。以上不含导入、部署、传输、后处理、制图和文档，未测部分不虚构总耗时。

共12条短轨迹：3条150步试跑及6条正式OFF/弱2 ms完成，3条强控制正式失败。停止原因是科学质量，不是耗尽10–15分钟预算；没有为追求成功延长。远端62个文件共7,928,266 bytes回传逐一SHA核对；源包、Supervisor配置/日志、全步轨迹、失败提案和输入身份保留。

四组真实PNG/PDF位于results/straight_tube_exp1/f10o/figures：Fig_F10O_1_contact_mechanism、Fig_F10O_2_particle_convergence、Fig_F10O_3_XY_navigation、Fig_F10O_4_admission_cost。脚本在scripts/straight_tube_exp1/f10o/figures.py；FIGURE_PROVENANCE.json逐输入和图件SHA、VISUAL_REVIEW.json实际检查记录一起归档。排版修订不改变计算数据。

科研进度文档只在验证完成后追加一段，备份为microbubble_aggregation.md.bak_F10O_20261013；before/after、diff、SHA和标题结构检查见progress_document/。

## 下一轮

先在现有低驱动可通过区域内登记少量中间幅值/相消控制，保留原刚度及重叠阈值，把真实泡处声压和每步接触准入作为控制约束；先跑2 ms及步长/缓存对照，再尝试X/Y独立短目标与小段反馈。形态目标应依据B2实际去质心张量，不能继续沿用M0的0.42%收益估计。若需要改变接触物理，必须单独登记依据和跨模型对照，不能作为本轮失败的补救性调参。不重写求解器、不加新CFD/材料模型、不开展AI训练。
"""
    (ROOT / "docs/FAST_PZT_F10O_DEVELOPMENT.md").write_text(text)
    (ROOT / "docs/FAST_PZT_F10O_B2_VALIDATION.md").write_text("""# F10O B2 validation

The engineering bridge and mathematical equations pass. Strong-control overlap admission fails at the 0.175 ms proposal for both time steps and both caches. All three valid prefixes and rejected proposals are retained; no stiffness or tolerance changes.

The adapter directly imports original F10H geometric_contact and F10G resistance. R3 GPU walls inspect both predictions and final proposals, and any correction rejects the interior prototype. Pressure/APECSS violations also fail closed. M0 uses two local lubrication iterations and PBD; B2 uses the same resistance model solved as a matrix and the geometric unilateral spring. Cross-model contrasts are not numerical error estimates.

Authoritative evidence: analysis/FINAL_GATE.json, COMPARISONS.json, MATCHED_0p1MS.json, RUN_ANALYSIS.json, each run's all-step trajectory.h5 and LAST_STEP.json. Topology reuses the historical persistence implementation and all accepted-step edges; sub-2ms prefixes cannot establish a 2ms sustained event.

Weak-control two-step/two-cache morphology is consistent within the tested discrete sensitivity but is much smaller than historical M0. This does not calibrate contact, near-field acoustic response or wall hydrodynamics. Full Chinese results and exact identities: FAST_PZT_F10O_DEVELOPMENT.md.
""")
    (ROOT / "docs/FAST_PZT_F10O_XY_NAVIGATION.md").write_text("""# F10O XY navigation admission

PREREGISTRATION.json locked X=(0.2,0,0) µm and Y=(0,0.2,0) µm, tolerance 0.075 µm, before the pilot. The historical F10N (3,2,1) µm target is unchanged.

Status: NOT_RUN_STRONG_B2_QUALITY_GATE_FAILED. New XY baseline, optimization and feedback were not launched. Numerical OFF/weak 2ms controls and the strong valid prefix are labelled as such in Figure 3. Neither reduced error over the truncated prefix nor the previous F10N COM displacement is claimed as new XY target success. Independent X/Y control remains untested. Symmetric near-zero Z COM is not a proof of universal Z uncontrollability.

5, 10 and 30ms windows are NOT_RUN. Next: first establish a lower-drive admissible 2ms control window with full fixed-ID sensitivity checks, then preregister independent X/Y policies and bounded feedback. No contact-budget relaxation.
""")
    write(
        BASE / "analysis/DELIVERY_IDENTITIES.json",
        dict(
            baseline=g.get("baseline_commit", "e2dfb4f41cea33082583bbf8424a963ee194e85c"),
            cleanup="19c11e885",
            bridge="3f5a117c1",
            scientific_runs=commits,
            report_sha256=sha(ROOT / "docs/FAST_PZT_F10O_DEVELOPMENT.md"),
        ),
    )


if __name__ == "__main__":
    main()
