# V3 chunk fragmentation mechanism audit

Task: `V3_CHUNK_FRAGMENTATION_MECHANISM_V1`, run `20261015T174630Z`.

The actual scientific conclusions and run status are in:

- `verification/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z/FINAL_STATUS.json`
- `verification/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z/REVIEW_V3_CHUNK_FRAGMENTATION_V1.md`
- `results/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z/index.html`

The historical V3 run is read-only. Its full manifest, native binary, original identities, accepted node stream, lossless trial archives and sparse damage/contact events are required on the server. The GitHub snapshot intentionally omits some of these inputs. Do not fabricate missing every-step node histories or interpolate rejected trials into accepted trajectories.

The server environment remains `/workspace/y_clot/env.sh`, project root `/workspace/y_clot/project`; `Y_CLOT_PROJECT_ROOT` overrides the root for offline diagnostics. The native library is the unchanged frozen V3 `liblammps.so`, SHA256 `da705cdd2358eab259d9f99ab402e743f6a64c829959eb094c148f8492552ae7`. Native runs additionally put its build directory first in `LD_LIBRARY_PATH`.

Processing order: `prepare.py`, `topology_audit.py`, `recover_states.py`, `spatial_audit.py`, `mixed_mode.py`, `group_audit.py`, `augment_audits.py`, `benchmark.py`, `benchmark_analysis.py`, `figures.py`, `finalize.py`, actual WSL document insertion, and `seal.py`. The shared `graph_ops.py` uses original live bulk/interface bonds and never counts contact as a healed material connection. `bench_native.py` loads the exact frozen bridge classes and replaces only the synthetic reference geometry and prescribed boundary motion.

Do not rerun `prepare.py` or the frozen benchmark into an existing result directory. They deliberately refuse to overwrite frozen designs. A fresh authorized replication requires its own run identifier and explicit baseline paths. The corrected benchmark also verifies the prior accepted prefix before continuation; its original failed interface setup and rejected grouped replay remain archived. A fresh independent replay should use the delivered design cards, not silently remove those checks.

`ACTIVE_DEVELOPMENT_ENTRYPOINTS.md` at the project root is itself included in the protected historical V3 SHA manifest. This task therefore preserves it byte-for-byte and provides this new entrypoint instead. No Git push is part of this development task.

Scientific boundaries: kinematic shear is a diagnostic, single-bond virtual removal uses the complete endpoint-star minimized energy, graph weight `1-d` is not certified load capacity, and no reported U/D/R/virtual energy is calibrated material Gc. Synthetic specimen results and actual V3 fragmentation remain separate.
