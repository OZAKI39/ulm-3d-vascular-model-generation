# CFD validated evidence index

This is the compact human-readable index. Machine-readable paths, statuses,
and hashes are in `cfd_validated_evidence_index.json` in the slimming QC run.

| Topic | Status | Canonical implementation/test | Canonical evidence | Retired alternative |
|---|---|---|---|---|
| Tau=1 scaling and pressure | Validated | `validated_contract.py`; `test_validated_contract.py` | Base final QC | fixed 23622 Pa and mis-scaled runners |
| Geometry and continuous q | PASS | `dimensionless_geometry_kernel.py`; scale tests | final Base geometry validation | fallback-q repair runners |
| Adaptive inlet population | PASS, 375 active | `adaptive_target_population.py`; its test | root-cause and fix JSON | 376-element denominator/debug runners |
| Full-timestep conservation | PASS | `full_timestep_mass_referee.py`; Full V2 tests | Base Full V2 and final contract JSON | boundary-only final referee |
| Physical port planes | PASS | `physical_port_flux.py`; geometry contract test | V3 plane JSON, SHA `ffaa49...f005` | standardized-plane one-shot writer |
| Physical volume flux | Validated | `physical_port_flux.py`; analytic MLS tests | Base/Coarse final QC | CELL_CUBE V1 |
| Base steady | `CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS` | `steady_state.py`; evidence loader tests | Base final QC and restart `ffcd98...c39f` | Base long-run orchestration |
| Coarse sensitivity | PASS two-grid only | optional `tau1_grid_convergence.py` | Coarse final and two-grid JSON | C/B/F solver runner |
| Fine | `NOT_COMPLETED_RESOURCE_BUDGET` | adaptive regression/evidence loader | root cause, fix, 5,000 safety, termination | Fine monitors and million-step continuation |

Key files are loaded by `utils/cfd_flow/validated_evidence.py`. Large restart
and mesh payloads are addressed by path, size, and SHA256 and should not be
opened as text. Formal three-grid GCI remains incomplete.
