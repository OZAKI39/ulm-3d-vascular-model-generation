# F10Q active code map

Complete F10P baseline: `1c785f4803ca7dc0ea2fe91756b02487c7c75a59`, clean and independently verified. Separate F10Q worktree.

- Four-PZT complex field / analytic cubic gradient: `acoustics/water_disc_proxy.py`, `fast_sim/cubic_field.py`.
- Frozen F10G pair resistance and geometry: `fast_sim/f10g_contact.py`.
- Protected F10P local QP and adapter: `fast_sim/hard_sphere.py`, `water_hard_sphere.py`.
- Opt-in F10Q sequential solver / adapter: `fast_sim/sequential_hard_sphere.py`, `water_sequential.py`.
- Original APECSS / F2 / GPU walls: `water_b2.py`, `water_control_engine.py`, `nearfield.py`, `models.py`, `pdms_cavity_engine.py`.
- B2 / M0 references: unchanged `f10h_contact.py` / `water_control_engine.py`.
- Fixed-ID group metrics and layered gates: `scripts/straight_tube_exp1/f10o/analyze.py`, `f10p/analyze.py`.
- F10Q entrypoints: `scripts/straight_tube_exp1/f10q/`; protocol: `configs/straight_tube_exp1/f10q/`.
- Exact prior inputs and failed proposals: `results/straight_tube_exp1/{f10n,f10o,f10p}`; never overwritten.

## Bounded cleanup

本轮无新增可证明安全删除的对象。Actual reference search and protected hashes are stored in `results/straight_tube_exp1/f10q/cleanup/`. Conditional historical CLI entrypoints remain reproduction dependencies. No numerical source was changed in the separate cleanup audit commit. F1–F10P results, failures, source snapshots and all historical implementations are retained.
