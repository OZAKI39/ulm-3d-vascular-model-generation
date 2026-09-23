# Particle-8：连续混合入口可视化与回放

本阶段只增加回放、图像与动画导出，不修改 Frozen FEM 或 P0–P7 科学模型。用户已在本阶段指令中接受 P6.5/P7 为基础；仓库原历史审核快照保持不变。

真实证据与合成证据分别保存为 `real_frozen` / `synthetic_control` 场景。真实单微泡沿 P7 的 5 个 P6.5 接受状态显示插值；真实混合 smoke 中唯一已接纳的 RBC 维持静态障碍，不从诊断速度伪造运动。排队粒子没有管腔位置。RBC 连续 admission/passage 仍未建立。

合成动画复用 P7 实际 LAMMPS 二进制重启验证的 2358 个事件，而不是重跑 110 万粒子长账本。控制区长 10 nm、宽 200 µm，速度 25 µm/s。有限形状跨越开口，中心在区内驻留 0.4 ms。左图流向坐标拉伸、ID 符号不代表物理形状；右图显示原尺寸和原姿态的真实 xy 投影。三个慢放物理窗口与末端记账之间有明确剪切，不改变任何事件时刻。

固定前提：H_D=0.45 为 feed volume fraction；C_MB=8.5e12 m⁻³ 为名义估算；确定性累计通量调度；各向同性 RBC 取向为模型假设。tube Hct 不被强制为 45%。生产 timestep/neighbor 参数及非球形润滑未冻结，P6.5 事件时间收敛未建立。没有 full suspension、完整 PK、CFD 或 Particle-7.5。

## 复现

在项目根目录运行，使用现有 `.venv`。额外依赖只提供本地 FFmpeg 编码/解码：

```bash
.venv/bin/pip install -r particle_3d/requirements-particle8.txt
.venv/bin/python particle_3d/scripts/export_particle8_replay.py
.venv/bin/python particle_3d/scripts/render_particle8.py
.venv/bin/python particle_3d/scripts/verify_particle8.py
.venv/bin/python particle_3d/scripts/finalize_particle8.py --decode --final-gate
```

`export_particle8_replay.py` 首先核验原 P7 源码及数据哈希，仅从上游缓存读取，不调用求解器或随机采样器。导出 schema `PARTICLE8_REPLAY_V1`、原样本、geometry、稳定 ID、全部生命周期事件与轨迹。SI 数据为权威值，出生时刻定义为实际接纳，预约时刻另存。事件状态保存 `scheduled/pending/admitted/active/exited/deleted`；帧状态右连续，接纳后活动，退出后同一事务删除，不在未来帧保留幽灵粒子。

只重画某个动画：

```bash
.venv/bin/python particle_3d/scripts/render_particle8.py --animations 3 --animations-only
.venv/bin/python particle_3d/scripts/render_particle8.py --figures-only
```

可以为 export/render 指定 `--output /tmp/particle8-replay`，生成独立副本。MP4 为 1440×900 / 15 fps，编码 H.264、yuv420p、faststart。GIF 是 3 fps 缩小预览。帧率和播放窗口只是显示设置，不是生产物理 timestep。Flux 动画使用显示时钟，候选点没有被赋予虚假的物理出生时刻。

`finalize --final-gate` 要求所有测试通过、视频完整解码、代理视觉检查凭据与产物哈希一致。重新渲染若改变像素或视频，需要重新目检并更新 `data/visual_inspection.json`，不能沿用失效凭据。此文件中的代理 PASS 不等于用户验收；用户状态单独保留 `PENDING_USER_REVIEW`。

## 入口与产物

打开 `reports/particle8/index.html` 可浏览所有 MP4、GIF 和 PNG；Markdown 报告为 `PARTICLE8_REVIEW.md`，机器记录为 `PARTICLE8_VALIDATION.json`。`OUTPUT_FILES.json` 列出全部非日志输出的相对路径、大小和 SHA256；`logs/` 保存命令、软件测试和 JUnit。每张图有来源映射，每部动画有完整状态的逐帧 manifest，视频解码联系表用于审核。
