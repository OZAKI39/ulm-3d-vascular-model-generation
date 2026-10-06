# Native P2/P1 staged validation

Start with [the Chinese review](NATIVE_P2P1_DEVELOPMENT_REVIEW_ZH.md) or [figures](OPEN_RESULTS.html).

Final state: `P2P1_NOT_EVALUABLE_DUE_TO_IMPLEMENTATION_BLOCKER`. Affine representation and reference sampler pass; native curved mixed geometry does not pass the pre-CFD gate. No production CFD or trajectories were run.

These scripts are isolated validation tools. They use the project root inferred from this directory; source/build inspection additionally uses the retained production build identities. The CPU sampler is not installed into the particle package.

Completed supervised stages, in dependency order:

1. `preflight.py`: verify retained source/binary/PETSc and protected inputs, record actual resources. `runtime_identity.py` also rechecks CUDA driver/runtime/cuSPARSE/PETSc library hashes without GPU computation.
2. `native_build.py`: link a unit harness against production objects, save native basis samples.
3. `path_audit.py`: bind the inspected call chain to native source hashes.
4. `representation.py`: B1/B2 analytic fixtures and real B3 WALL-edge topology.
5. `sampler_tests.py`: validation-only NumPy sampler versus independent analytic and native basis values.
6. `poiseuille_representation.py`: deterministic affine elevation of frozen V2, all original probes and random points.
7. `curved_geometry_audit.py`: two native curved single-cell assembly probes, with affine controls. This triggers the stop gate.
8. `figures.py`, `final_tests.py`: saved-data plots and regression checks.
9. `finalize.py`: short receipt/identity/gate packaging after all supervisors exit.

For an explicitly authorized reproduction, read the server operating guide and inspect resources first. Each scientific stage uses the existing project supervisor through `python supervise.py <stage>`. The non-CFD cumulative budget is 5400 s, including failed attempts. The wrapper does not enqueue a next stage. A failed science gate is a valid process exit and must be read from its JSON before starting a dependent stage. There is no CFD launcher in this directory.

`runtime/PIPE_AFFINE_ELEVATION.npz` is retained on the compute server and reconstructible from the frozen mesh without a mesher; its SHA is in the identity. Runtime harness binaries, build products, and detailed monitoring logs remain outside Git. Compact build/object identities, native basis data, geometry samples, figure data and all stage receipts are preserved in tracked reports/NPZ. `OUTPUT_SHA256.json` covers the final delivered artifacts and excludes itself, runtime, monitoring and TASK state.

Tests certify engineering behavior and detect the curved assembly inconsistency. They do not certify CFD or trajectory accuracy.
