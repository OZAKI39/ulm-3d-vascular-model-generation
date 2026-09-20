# Frozen 3D FEM / SimVascular baseline

FEM DEVELOPMENT: **FROZEN after Stage SV1.3Q**。本目录保存粒子开发需要的源码、accepted mesh/flow/boundaries、配置、证据与校验；尚未实现 Particle-0。

从 [PARTICLE_HANDOFF.md](PARTICLE_HANDOFF.md) 开始，再读 [FROZEN_FEM_BASELINE.md](FROZEN_FEM_BASELINE.md) 和 [PARTICLE_RESEARCH_ROADMAP.md](PARTICLE_RESEARCH_ROADMAP.md)。同步详情见 [SYNC_REPORT_FEM_PARTICLE_HANDOFF.md](SYNC_REPORT_FEM_PARTICLE_HANDOFF.md)，历史及可移植性见 [PORTABILITY_NOTES.md](PORTABILITY_NOTES.md)。

新环境：Python≥3.11，`python -m pip install -e .`。只读验证：

```bash
python -B scripts/fem_freeze_sync/validate_frozen.py
python -B scripts/fem_freeze_sync/verify_manifest.py
```

科学数值源于实际 Stage Q artifact，不源于聊天手抄。Mesh convergence / time-step sensitivity 未做且不计划在粒子开发前补做；CPU/GPU field equivalence deferred，均为用户接受的项目决定，不是 PASS。禁止把本分支当作开始进一步 FEM 调优的授权。
