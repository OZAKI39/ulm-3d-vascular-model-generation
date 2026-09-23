# Particle-8 完整三维可视化增强：中文审核报告

本阶段已完成显示层代码、永久测试、五组 MP4/GIF、关键帧、storyboard 和机器记录。审核结论仅适用于已有证据的回放展示，真实 RBC 连续通行仍未建立。

```
AUTOMATED_CHECKS = PASS
VISUAL_OUTPUTS = PASS
SCIENTIFIC_LABEL = PASS
REPLAY_FIDELITY = PASS
STAGE_RESULT = PASS_WITH_LIMITATIONS
```

代码提交：`d2c7b82e65467c2c44fd1cb62c95fd72d1fee9fd`。新输出：`/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/outputs/particle8_full3d`。既有 1677 个受锁文件 SHA256 保持一致。

本轮实际运行 68 项测试（新增 full3d 32 项，原 P8 36 项），全部通过；895 个 MP4 帧和 179 个 GIF 帧均成功解码。先前 P0–P8 的 641 项是继承证据，本轮未重跑整套物理测试。

初次检查发现离散胶囊包围盒中心不是物理中心，测试已改为解析形状表面约束；同时改进 MB 显示网格的双精度坐标构造。抽查视频发现 A 出生前误用合成 section 说明，已修正分支、补回归测试并重新导出。早期 `logs/core_tests.*` 保留为开发记录，最终以 `logs/full3d_tests.*`、`logs/upstream_particle8.*` 和 `data/test_run.json` 为准。

自动发现并复用原 P8 四个 scene、事件/生命周期账本、稳定 ID、原始几何、姿态、物理时钟和重启缓存；复用 P8 snapshot / interpolate 与 figure 配色。Frozen 模型从官方 boundary manifest 读取，主视图包含全部 45,221 个 wall 三角形以及 INLET、OUTLET_01/02/03，无抽稀和局部截取。

渲染：PyVista/VTK 离屏、正交投影；MP4 H.264/yuv420p、1600×1000、15 fps，GIF 960×600、每五帧采样、约 3 fps。GIF 时间量化用于预览，以画面 physical clock 和 MP4 帧记录为准。

相机：主视图方位角按视频帧线性从 35° 增到 115°，俯仰恒定 28°；速度约 6.2–7.5°/视频秒。真实入口副图比主视图多 90°。固定 A 主视图为 35°、副图 125°，两版本物理状态散列完全相同。相机、原点平移和 µm 单位换算均仅用于显示，不旋转物理坐标。

插值：复用 P8 线性位置与 SLERP 姿态，仅用于显示；不外推额外物理解、不改变 birth/admit/delete 时序。尾迹为保存折线加近期窗口和当前时刻的裁切插值端点，无样条平滑。默认真实 0.35 ms、合成 0.25 ms，只有 active 粒子有尾迹，删除时移除；静止 RBC 不制造非零轨迹。支持 `--tail-mode full` 另行导出，但本交付审核针对 recent。

真实粒子保持原始尺寸与姿态；capsule 使用保存的 capsule_axis。A 放大图以原半径淡色轮廓、中心标记显示 MB，让不足一个半径的运动可见。合成主图 z×20,000 且采用中心 ID glyph，均显式标记 display-only；C 副图显示未拉伸原形，D 副图同尺度显示 restarted 对照。合成中心轨迹不代表有限形状完全穿过 section。

## 逐动画审查

### A_real_single_mb_orbit

- 数据源：`particle_3d/reports/particle8/data/real_single_mb_scene.json`；真实 Frozen / 单 MB；保存状态来自 P6.5，经 P7/P8 回放。
- 物理范围：仅入口附近 1 ms、5 个保存状态，位移约 0.440609 µm。主视图完整血管，副图保持原尺寸放大观察。
- 显示插值：是，仅保存状态之间。
- 相机：orbit；尾迹：recent，0.35 ms。
- 限制：不能解释为贯穿完整血管；与 B 是独立 fixture，不能拼成真实混合输运。出生前无粒子。
- 输出：[MP4](mp4/particle8_full3d_A_real_single_mb_orbit.mp4) · [GIF](gif/particle8_full3d_A_real_single_mb_orbit.gif) · [storyboard](storyboard/particle8_full3d_A_real_single_mb_orbit_storyboard.png)。
- 关键帧：[首帧](keyframes/particle8_full3d_A_real_single_mb_orbit_decoded_0000.png) · [中帧](keyframes/particle8_full3d_A_real_single_mb_orbit_decoded_0090.png) · [末帧](keyframes/particle8_full3d_A_real_single_mb_orbit_decoded_0179.png)。

### A_real_single_mb_fixed

- 数据源：`particle_3d/reports/particle8/data/real_single_mb_scene.json`；同旋转 A，物理帧逐项一致，固定相机对照。
- 物理范围：同一 1 ms 单 MB 保存轨迹；固定视角便于区分短距离位移与相机引起的投影变化。
- 显示插值：是，仅保存状态之间。
- 相机：fixed；尾迹：recent，0.35 ms。
- 限制：固定与旋转版都只回放已有状态，不证明新时间步收敛或全血管通行。
- 输出：[MP4](mp4/particle8_full3d_A_real_single_mb_fixed.mp4) · [GIF](gif/particle8_full3d_A_real_single_mb_fixed.gif) · [storyboard](storyboard/particle8_full3d_A_real_single_mb_fixed_storyboard.png)。
- 关键帧：[首帧](keyframes/particle8_full3d_A_real_single_mb_fixed_decoded_0000.png) · [中帧](keyframes/particle8_full3d_A_real_single_mb_fixed_decoded_0090.png) · [末帧](keyframes/particle8_full3d_A_real_single_mb_fixed_decoded_0179.png)。

### B_real_mixed_smoke_orbit

- 数据源：`particle_3d/reports/particle8/data/real_mixed_scene.json`；真实 Frozen / mixed inlet smoke；仅绘制实际 admitted 且 active 的粒子。
- 物理范围：0–42.986308099 s。结尾 RBC scheduled=1110、admitted=active=1、pending=1109；MB scheduled=pending=1、admitted=0。
- 显示插值：是，仅保存状态之间。
- 相机：orbit；尾迹：recent，0.35 ms。
- 限制：唯一 admitted RBC 是保存的静止 capsule surrogate。尾迹长度为零，不根据诊断速度制造运动。pending 从不画入血管。
- 输出：[MP4](mp4/particle8_full3d_B_real_mixed_smoke_orbit.mp4) · [GIF](gif/particle8_full3d_B_real_mixed_smoke_orbit.gif) · [storyboard](storyboard/particle8_full3d_B_real_mixed_smoke_orbit_storyboard.png)。
- 关键帧：[首帧](keyframes/particle8_full3d_B_real_mixed_smoke_orbit_decoded_0000.png) · [中帧](keyframes/particle8_full3d_B_real_mixed_smoke_orbit_decoded_0090.png) · [末帧](keyframes/particle8_full3d_B_real_mixed_smoke_orbit_decoded_0179.png)。

### C_synthetic_lifecycle_orbit

- 数据源：`particle_3d/reports/particle8/data/synthetic_scene.json`；SYNTHETIC CONTROL / NOT REAL FROZEN LUMEN；P8 synthetic scene。
- 物理范围：宽 200 µm、长 10 nm 的开放 section。三个连续观察窗口与最终记账帧之间有明确 CUT 标签。
- 显示插值：是，仅保存状态之间。
- 相机：orbit；尾迹：recent，0.25 ms。
- 限制：最终 2356 RBC 与 2 MB admitted；2352 RBC 与 2 MB exited/deleted；4 RBC active。pending=0 是该控制案例的实际结果，未虚构排队。有限形状跨越开放端盖，不是生理通过证据。
- 输出：[MP4](mp4/particle8_full3d_C_synthetic_lifecycle_orbit.mp4) · [GIF](gif/particle8_full3d_C_synthetic_lifecycle_orbit.gif) · [storyboard](storyboard/particle8_full3d_C_synthetic_lifecycle_orbit_storyboard.png)。
- 关键帧：[首帧](keyframes/particle8_full3d_C_synthetic_lifecycle_orbit_decoded_0000.png) · [中帧](keyframes/particle8_full3d_C_synthetic_lifecycle_orbit_decoded_0097.png) · [末帧](keyframes/particle8_full3d_C_synthetic_lifecycle_orbit_decoded_0194.png)。

### D_restart_continuity_orbit

- 数据源：`particle_3d/reports/particle8/data/synthetic_scene.json`；SYNTHETIC CONTROL / RESTART AUDIT；continuous 与 restarted 两份 P8 scene。
- 物理范围：0.124–0.126 s，checkpoint=0.125 s，并切到 0.25 s 最终记账。两视图采用相同拉伸和相机。
- 显示插值：是，仅保存状态之间。
- 相机：orbit；尾迹：recent，0.25 ms。
- 限制：所有显示帧状态散列相同、mismatch=0。使用已接受的 P7 binary restart 证据缓存，本轮没有重跑物理重启或产生新解。
- 输出：[MP4](mp4/particle8_full3d_D_restart_continuity_orbit.mp4) · [GIF](gif/particle8_full3d_D_restart_continuity_orbit.gif) · [storyboard](storyboard/particle8_full3d_D_restart_continuity_orbit_storyboard.png)。
- 关键帧：[首帧](keyframes/particle8_full3d_D_restart_continuity_orbit_decoded_0000.png) · [中帧](keyframes/particle8_full3d_D_restart_continuity_orbit_decoded_0080.png) · [末帧](keyframes/particle8_full3d_D_restart_continuity_orbit_decoded_0159.png)。

## 保留的科学边界

Real RBC passage remains unresolved。Synthetic control 仅用于 bookkeeping / lifecycle；真实 smoke 的静止 RBC 不能解释为连续输运。Production timestep 与 neighbor settings 未冻结。近场仍为 sphere-normal regularization，非球形 lubrication、切向/旋转耦合与完整多体 mobility 未冻结。H_D=0.45 是 feed volume fraction，不是瞬时 tube hematocrit；tube Hct 仅为独立诊断。没有 Particle-7.5、CFD、新 RBC 变形物理、full suspension 或 PK。原始 P8 全部 carry-forward limitations 原样保存在验证 JSON。

## 审核证据与范围

所有 MP4/GIF 已逐帧解码；逐帧数值审计核对位置、q、几何散列、stable IDs、pending 排除与 restart 一致性。代理视觉检查查看每组实际视频解码的首/中/末帧 contact sheet、中帧原尺寸图片、总览和科学范围图；未声称人工逐帧看完整视频，也不代替用户验收。用户验收状态为 PENDING_USER_REVIEW。查看文件和视频 SHA256、每组观察结论见 `data/visual_inspection.json`。

- [本地浏览与 storyboard 页面](index.html)
- [总览 PNG](storyboard/full3d_storyboard.png) · [科学范围总结图](storyboard/scientific_scope_summary.png)
- [机器验证](PARTICLE8_FULL3D_VALIDATION.json) · [测试执行记录](data/test_run.json) · [解码审计](data/media_audit.json)
- [逐帧 CSV](data/frame_audit.csv) · [active 原始坐标/姿态 CSV](data/active_particles.csv)
- [完整文件清单](OUTPUT_FILES.json) · [产物 SHA256](ARTIFACT_SHA256.json)
