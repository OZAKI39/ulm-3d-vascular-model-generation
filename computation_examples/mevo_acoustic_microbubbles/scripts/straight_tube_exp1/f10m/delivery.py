"""Current evidence report and explicit engineering/science admission split."""

from .common import BASE, ROOT, read, write, identity


def main():
    g = read(BASE / "analysis/GEOMETRY.json")
    r = read(BASE / "analysis/RESOURCE_PREFLIGHT.json")
    a = read(BASE / "adapter_validation/VALIDATION.json")
    b = read(BASE / "source_bridge/GATE.json")
    d = read(BASE / "dynamics/SUMMARY.json")
    cpu = read(BASE / "analysis/CPU_CHECKS.json")
    transfer = read(BASE / "analysis/TRANSFER_VERIFICATION.json")
    config = read(BASE / "CONFIG_SNAPSHOT.json")
    assert all(q["status"] == "PASS" for q in [g, a, b, cpu, transfer])
    assert d["status"] == "COMPLETE" and all(d["checks"].values())
    assert g["end_conditions"] == ["open_fail_stop", "closed"] and g["top_port"]["diameter_m"] == 0.0015
    review = read(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS" and len(review["figures"]) == 4
    gate = dict(
        engineering_status="PASS",
        full_apparatus_acoustics="NOT_RUN",
        three_basis_count=0,
        ON="NOT_RUN",
        OFF="COMPLETE_1_MS_MAIN_BORE_ONLY",
        junction_particle_motion="NOT_ADMITTED_FAIL_STOP",
        geometry="CONFIRMED_BLOCK_AND_MAIN_BORE_WITH_TOP_PORT; port reference/centering/blind depth provisional",
        numerical_checks=dict(
            CPU=len(cpu["checks"]),
            planar_GPU=len(a["checks"]),
            CW_bridge=b["status"],
            OFF_GPU=len(d["checks"]),
        ),
        **config["protocol"]["calibration"],
        identity=identity(),
    )
    write(BASE / "analysis/FINAL_GATE.json", gate)
    total = a["actual_wall_s"] + b["wall_s"] + d["whole_case_wall_s"]
    report = f"""# F10M：PDMS居中通道、三枚7 mm PZT与多介质接口

本轮工程验证通过；真实装置三路底稿 **NOT_RUN**，ON轨迹 **NOT_RUN**。这不是F10L均匀水场的重命名：已完成32.5×9×9 mm PDMS材料分区和新的多介质接口，以及实际GPU小型验证、真实主通道OFF轨迹。继续EXPLORATORY / MATERIAL_ACOUSTICS_UNCALIBRATED / SOURCE_PRESSURE_UNCALIBRATED / PHYSICS_UNCALIBRATED，training_ready=false、experiment_matched=false。

## 基线与输入身份

基于完整F10L {identity()["baseline_commit"]}，独立分支dev/fast-pzt-f10m-pdms-three-pzt-20261013。原F10L工作区和科学文件保留。当前工作区：{ROOT}。远端独立目录/workspace/mevo_acoustic_microbubbles/f10m_pdms_block，每次Supervisor执行来自.execution/<commit>不可变快照，执行SHA和源码哈希见execution及每项结果identity。未变更CUDA、Warp或项目软件环境，没有重跑旧CFD、旧声場或历史研究。

## 最新装置与本轮中途确认

PDMS边界X±16.25、Y/Z±4.5 mm，通道严格居中、沿X、内径1.6 mm、长度32.5 mm。没有独立薄管壁和外径。两侧PDMS路径3.7 mm。SMD07T05R412WL三枚有效陶瓷面均为7 mm：

| ID | 中心/mm | 法向 |
|---|---|---|
| side_y_minus | (0,−4.5,0) | (0,+1,0) |
| side_y_plus | (0,+4.5,0) | (0,−1,0) |
| end_x_plus | (+16.25,0,0) | (−1,0,0) |

无Z侧源。居中端片覆盖主通道开口100%，覆盖面积2.010619 mm²；不能用透明声源隐藏这一事实。用户随后确认：顶部开口直径1.5 mm，向下连通中心通道，注入后保持开放，X+端始终封闭。已更新为顶部圆柱孔与主通道的并集，端部粒子边界为封闭；修正前配置、轨迹和图件完整保留在history/pre_top_port。按“端PZT发声面到孔中心”暂定5.4 mm参考，孔轴x=10.85 mm，y=0暂按居中；没有擅自增加孔底盲段。孔的具体测距参考、横向中心与钻孔深度仍待核对；X−端暂维持原开放计算假设。

完整材料图包含顶部孔。逐泡模型本轮只准入远离交汇口的主通道：进入孔轴±(0.75 mm+泡半径)轴向邻域会报错停止，保留有效前缀，不把真实孔口当作圆柱壁反射。不是已经实现复杂交汇口的接触模型。

## 文献材料与声学接口

材料参考为Xu等2020论文（https://doi.org/10.1103/PhysRevApplied.13.054069）表II/III中的SYLGARD184、10:1、85°C组：ρ=1031.4 kg/m³，5 MHz纵波相速1029 m/s，α=214.8 Np/m，幂指数1.494。测量25°C、样品静置10天；固化时长未知。论文原文表格和公式已实际核对，配置不冒充用户品牌/固化实测。覆盖3–7 MHz；2.5 MHz不作静默材料外推。内外液体的1500 m/s、1000 kg/m³及零吸收均为明确水演示假设。

新适配器使用与网格同形状的sound_speed、density、alpha_coeff数组或带SHA的HDF5。仅允许一个全局alpha_power；转换单位为dB/(MHz^y cm)，零频声速按论文色散关系转换，避免相速与求解器色散重复计入。已在安装的k-wave-python0.6.2/Python/CuPy实际运行验证。仍为纵波等效流体，没有剪切波、完整弹性、PZT实体反射/背衬/胶层或电压声压传递标定。新CW步长按最大声速，传播等待按最小声速；保留Re[P exp(+iωt)]及F6一致梯度。

## 资源准入与实际GPU验证

完整装置候选频率3 MHz，网格间距55 µm，最短波长约{r["minimum_ppw"]:.3f}格；物理网格{r["shape"]}，外置16格PML后{r["expanded_shape"]}。BLI容差0.01对应支持半宽ceil(1/(πtol))=32格，计划边距36格；全装置实际BLI数组未生成，不能声称它已实际通过求解。按F10L实测15100 MiB缩放并加多介质余量，估计{r["estimated_peak_GiB"]:.2f} GiB、三路约{r["estimated_three_basis_seconds"]:.0f} s，超过20 GiB/600 s本轮预算和24 GiB设备，因此先停止完整场，不用粗网格强行运行。估计不是实测F10M装置显存；没有7 MHz或1–7 MHz扫描。

实际完成六个小型3D平面传播对照（横向周期、X外置PML），5 MHz：标量水/常数材料图、水–PDMS无损界面、PDMS无损/吸收、空间吸收界面。耗时{a["actual_wall_s"]:.3f} s，采样整卡峰值{a["sampled_device_peak_MiB"]:.0f} MiB。常数图与标量相对差{a["scalar_map_relative_error"]:.3g}；理论压力反射系数{a["reflection_expected_pressure_coefficient"]:.6f}、实际反射幅值{a["reflection_measured_magnitude"]:.6f}，幅值绝对差{abs(abs(a["reflection_expected_pressure_coefficient"]) - a["reflection_measured_magnitude"]):.6f}（约8.3%相对差，非界面收敛证明）；吸收幅比理论{a["attenuation_ratio_expected"]:.6f}、实际{a["attenuation_ratio_measured"]:.6f}；PDMS传播相位误差{a["pdms_phase_error_rad"]:.6f} rad。门槛在运行前登记，未放宽。

另一个微型三圆片CW接口检查实际耗时{b["wall_s"]:.3f} s、峰值{b["metadata"]["sampled_device_peak_VRAM_MiB"]:.0f} MiB；实际BLI边界净余19格，F6梯度相对差{b["F6_consistent_gradient_relative_error"]:.3e}。这是0.2 mm测试圆片，不是7 mm真实装置的缩尺预测，不可拿其208 Pa等局部数值当装置声压。三个幅相输入独立，使用空间材料图。初次微型检查误用了0.001 BLI容差，对应319格支持，边界预检拒绝且没有GPU求解；改回已登记的0.01后通过。失败目录/日志保留，没有放宽稳态、梯度或接触门槛。

## 真实OFF轨迹

最新几何下32泡、R=1.78 µm、固定ID，1 ms，dt=1 µs，101个保存帧，背景速度与梯度严格零，无声场加载。原2.5 MHz APECSS响应SHA为{config["protocol"]["bubbles"]["response_sha256"]}，仅用于合法历史参数OFF演示。3或5 MHz没有匹配表，不制造ON。主通道中心的OFF全部位置逐位不变，Rg={d["Rg_start_end_m"][0] * 1e6:.6f} µm、分量32，接受重叠和墙修正均0。积分及输出约{d["wall_s"]:.6f} s，含准备/检查单案例{d["whole_case_wall_s"]:.6f} s。新X+有限半径封闭约束、X−开放越界停止、孔口禁入保护和4096点严格零流都实际GPU检查通过。总计{len(cpu["checks"])}项CPU检查、{len(a["checks"])}项平面GPU检查、CW接口及{len(d["checks"])}项逐泡GPU检查通过。

当前接受的数值验证/轨迹分组耗时合计约{total:.2f} s（不含环境导入、失败预检、修正前OFF与本地绘图），不是整个研发耗时。无网格收敛、无粒子时间步扫描、无局部物理校准、无导航或形态控制成功声明。

## 图件、回传和复现

figures/下四组Fig_F10M_1_PDMS_apparatus、Fig_F10M_2_material_regions、Fig_F10M_3_acoustic_admission、Fig_F10M_4_no_flow_bubbles各有PNG/PDF，均已实际打开检查。图1/2使用最新顶部孔并显示X+封闭；图3明确完整装置NOT_RUN，只呈现实际小算例及准入；图4是真实OFF逐泡与质量统计。FIGURE_PROVENANCE与VISUAL_REVIEW保存来源、绘图源码、执行SHA、文件SHA。重要远程文件{transfer["files"]}项、{transfer["bytes"] / 1e6:.3f} MB已回传且SHA全部匹配；所有历史失败及修正前结果保留。

```bash
# WSL：在本项目目录，使用_shared/python_environment/bin/python
export PYTHONPATH=src:.
python -m scripts.straight_tube_exp1.f10m.verify
# 新的独立结果目录才可执行preflight、checks和GPU任务；正式执行前提交源码。
# remote模块将不可变快照部署并登记Supervisor，例：
# python -m scripts.straight_tube_exp1.f10m.remote validation --label NEW_UNIQUE_LABEL --timeout 300
# python -m scripts.straight_tube_exp1.f10m.remote source_bridge --label NEW_UNIQUE_LABEL --timeout 180
# python -m scripts.straight_tube_exp1.f10m.remote dynamics --label NEW_UNIQUE_LABEL --timeout 120
```

## 下一步

先确认顶部孔测距参考/横向中心/盲深和X−边界、PDMS品牌/配比口径/固化/温度、内外介质、实际单频、电压相位与声压标定、微泡半径分布及浓度。明确实体端PZT与孔口介质的声学边界；在有效频段内准备对应APECSS。通过经过验证的省内存后端或足够显存资源完成装置场，不能直接拿当前小算例导航。随后预登记远离孔交汇区的可分辨短目标与OFF/固定/反馈对照，证明可达性，再研究形态指标。交汇区粒子运动须单独验证，不能保证任意目标或形状。
"""
    (ROOT / "docs/FAST_PZT_F10M_DEVELOPMENT.md").write_text(report)


if __name__ == "__main__":
    main()
