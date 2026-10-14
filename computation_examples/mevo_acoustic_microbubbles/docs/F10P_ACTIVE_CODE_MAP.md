# F10P active code map

Complete baseline: F10O `98cf17b8ecd03e831628ed7005a8e1c83d354a40`. Read this map and F10O_ACTIVE_CODE_MAP.md rather than historical stage trees.

| Purpose | Reused source |
|---|---|
| Four independent sources and exact cubic gradient | acoustics/water_disc_proxy.py, fast_sim/cubic_field.py |
| Strict APECSS and F2 primary/secondary forces | fast_sim/water_b2.py: fields/control, water_control_engine.py, models.py, nearfield.py |
| Original symmetric pair resistance | fast_sim/f10g_contact.py: resistance, contact_geometry |
| Protected B2 reference | fast_sim/f10h_contact.py: geometric_contact; water_b2.py |
| Protected M0/PBD | fast_sim/water_control_engine.py and inherited engines |
| R3 cavity and GPU wall checks | coupling/pdms_cavity.py, fast_sim/pdms_cavity_engine.py, water_b2.py: wall_check |
| Fixed-ID shape and topology analysis | scripts/straight_tube_exp1/f10o/analyze.py; fast_sim/morphodynamics.py |
| Exact inputs, historical failures | results/straight_tube_exp1/f10n and f10o; never overwritten |

Planned opt-in additions: fast_sim/hard_sphere.py and water_hard_sphere.py; scripts/straight_tube_exp1/f10p/{run,analyze,figures,verify}. These are new entrances, not changes to historical defaults.

## Bounded cleanup decision

本轮无新增安全删除对象。F10O already removed its sole proved unused f64 helper. Targeted inspection of the current adapter, force, resistance, contact, wall and analysis imports shows functioning callers or historical reproduction dependencies. Standalone CLI functions without inbound imports are not dead code. Historical references, failed proposals, fields, response tables, snapshots and documents remain intact. No broad historical source scan or speculative pruning is warranted.

The separate stage-zero commit records this map, reference searches, baseline identity and successful original F10O read-only verification. No numerical source changes belong in that commit.
