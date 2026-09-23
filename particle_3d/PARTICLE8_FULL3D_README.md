# Particle-8 完整三维回放增强

这是独立的显示层。直接读取已接受的 Particle-8 场景和官方 Frozen 边界，不修改求解器、事件时序、稳定 ID、原始几何或姿态。新输出位于 `particle_3d/outputs/particle8_full3d/`，旧 `reports/` 被写入保护和 SHA256 锁保护。

新增模块：`particle8_full3d_data.py` 负责发现输入、完整边界导出、物理时钟和尾迹；`particle8_full3d_visuals.py` 负责 PyVista 离屏渲染、正交相机与 MP4/GIF 编码。依赖沿用仓库 `.venv` 的 numpy、PyVista/VTK、Pillow、imageio-ffmpeg，无新增物理依赖。

在仓库根目录执行：

```bash
.venv/bin/python particle_3d/scripts/prepare_particle8_full3d.py
.venv/bin/python particle_3d/scripts/render_particle8_full3d.py
.venv/bin/python particle_3d/scripts/audit_particle8_full3d.py
.venv/bin/python particle_3d/scripts/verify_particle8_full3d.py
```

审核实际视频提取的首、中、末帧后，在新输出的 `data/visual_inspection.json` 记录逐案例观察及实际查看文件的 SHA256。此记录必须来自实际检查，不能由渲染成功自动推导。最后执行：

```bash
.venv/bin/python particle_3d/scripts/finalize_particle8_full3d.py --source-commit <本阶段代码提交>
```

最终脚本要求测试源文件、产物、解码报告、审核图片的散列保持一致，并核对指定代码提交。输出中文审核报告、机器记录、CSV、HTML 浏览页和完整文件清单；不自动宣称用户验收通过。

默认产生五组 MP4/GIF：A 真实单 MB 旋转、A 固定对照、B 真实混合入口 smoke、C 合成生命周期、D 合成 restart 一致性。全部 MP4 为 1600×1000、15 fps，GIF 为 960×600、约 3 fps。`index.html` 是本地浏览入口。

尾迹默认仅保留 active 粒子的最近历史：真实案例 0.35 ms、合成案例 0.25 ms。可用 `--tail-mode full` 输出带 `_full_tail` 后缀的完整历史版本，不覆盖 recent 版本；全尾迹仍会在删除事件时随粒子移除。可用 `--cases A_real_single_mb_orbit` 选择案例，`--preview` 输出中帧预览。默认审核、测试和报告针对五组 recent 产物，额外 full 产物需要单独审核。预览不能替代实际编码视频的审核。

真实主视图包含原始 45,221 个 wall 三角形、INLET 和全部 3 个 OUTLET，未裁剪或抽稀。相机绕固定物理坐标缓慢旋转，主视图 azimuth 35°→115°、elevation 28°；真实入口放大图 azimuth 比主视图多 90°。固定 A 保持 35°，物理状态序列与旋转 A 完全相同。

显示插值复用 P8 的线性位置和姿态 SLERP，尾迹仅连接保存状态及裁切端点，不产生新物理解。合成主视图将 z 放大 20,000 倍，仅绘制中心 ID glyph；C 的副图显示未拉伸的原始有限形状，D 副图使用相同拉伸展示重启对照。真实坐标仅平移固定原点、转换为 µm，无轴向拉伸。RBC 胶囊使用原始保存轴，不误用保留的 q 作为胶囊轴。

真实单 MB 仅有 1 ms、约 0.440609 µm 的保存运动，完整血管展示不等于全血管穿行。真实 mixed smoke 仅有一个静止 RBC 障碍；pending 只计数。合成 10 nm 开放 section 是生命周期记账示例，有限粒子跨越开放端盖，不能证明真实生理通行。真实 RBC passage、production timestep、production neighbor、非球形 lubrication 仍未建立或冻结；近场仍为 sphere-normal regularization。H_D=0.45 是 feed volume fraction，tube Hct 单独诊断。无 Particle-7.5、CFD、新 RBC 变形物理、full suspension 或 PK。
