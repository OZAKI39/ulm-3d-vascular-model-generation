# F10O active code map

Baseline: complete F10N e2dfb4f41cea33082583bbf8424a963ee194e85c. Use this index instead of rereading all historical stages.

| Function | Active entry / necessary input |
|---|---|
| Four-source water pressure and consistent gradient | src/mevo_mb/acoustics/water_disc_proxy.py: CachedWaterBasis; fast_sim/cubic_field.py |
| Full particle M0 | src/mevo_mb/fast_sim/water_control_engine.py: WaterControlSimulator |
| CPU B2 | src/mevo_mb/fast_sim/f10h_contact.py: geometric_contact; f10g_contact.py: resistance; F10K scaling.py uses the same function |
| Cavity walls | src/mevo_mb/coupling/pdms_cavity.py; fast_sim/pdms_cavity_engine.py; sealed X and exposed union, no Y-vessel substitute |
| APECSS | src/mevo_mb/fast_sim/models.py: ResponseSection; nearfield.py: F2Physics/pair_reference |
| Existing trajectories | scripts/straight_tube_exp1/f10n/campaign.py; results/straight_tube_exp1/f10n/campaign_v1/runs |
| Existing quality admission | scripts/straight_tube_exp1/f10n/analyze.py and verify.py; fast_sim/f6_engine.py; f10h_diagnostic.py |
| Existing figures | scripts/straight_tube_exp1/f10n/figures.py |
| Required configurations | configs/straight_tube_exp1/f10n/{geometry,protocol}.yaml |
| Required data | results/straight_tube_exp1/f10n/{initial_states.npz,PREREGISTRATION.json,inputs/response.h5,cache/coarse.h5,cache/fine.h5}; full-precision runs/*/INPUT.json |

F10O entrypoints added after cleanup commit 19c11e885 and the passing F10N verifier:

- Adapter: src/mevo_mb/fast_sim/water_b2.py (calls unchanged geometric_contact and resistance).
- Run / deploy / archive: scripts.straight_tube_exp1.f10o.{run,remote,archive}.
- Analysis / quality / figures: scripts.straight_tube_exp1.f10o.{analyze,gate,figures}.
- Frozen protocol and targets: configs/straight_tube_exp1/f10o/protocol.json; results/straight_tube_exp1/f10o/PREREGISTRATION.json.
- Actual attempts: results/straight_tube_exp1/f10o/{pilot_v1,validation_v1}/.

Historical Y-vessel, three-source, old GPU B2 and retired experiment runners remain reproduction dependencies, not new control entrypoints.
