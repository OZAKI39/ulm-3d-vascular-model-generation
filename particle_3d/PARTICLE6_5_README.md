# Particle-6.5

Explicit sphere normal near-field V1. P0–P6 default numerical sources and P5 historical evidence stay immutable. No CFD or Particle-7.

The exact contract is [NEAR_FIELD_REGULARIZATION_V1.json](contracts/NEAR_FIELD_REGULARIZATION_V1.json). Raw geometry `h_geom`, coefficient denominator `h_eff=max(h_geom,h_lower)` and constraint gap `g_NF=h_geom-h_lower` are distinct. Original P5 leading coefficients and resistance metric are reused; original P3 physical-time subdivision certifies every accepted interval. Query envelopes may include h_lower for conservative candidate coverage; actual particle supports, radii, triangle meshes and closest points remain unchanged. A continuum handoff is not molecular solid contact.

Run from repository root with the P6 project environment:

```bash
.venv/bin/python particle_3d/scripts/run_particle65_validation.py
.venv/bin/python particle_3d/scripts/generate_particle65_report.py
.venv/bin/python particle_3d/scripts/finalize_particle65_report.py
.venv/bin/python -B -m pytest -q -p no:cacheprovider particle_3d/tests/particle6_5
```

The validation driver defaults to recomputing every trajectory; `--reuse-trajectories` is an explicit evidence-recovery option, not the numerical verification path. Figure generation reads saved evidence only; `--figures-only --output-dir PATH` performs isolated regeneration. Final report PASS requires successful final regression logs, not a hardcoded verdict.

V1 motion: `Particle65Stepper` in `particle65_motion.py`. A separate `particle65_checkpoint.py` stores the V1 contract, floor, original bridge schema, physical time, IDs, provenance and file hashes around an actual P6 LAMMPS binary restart. It rejects contract/provenance mismatches. P6 restart APIs are not changed.

The original P5 two-MB initial gap is already below 2 nm and is explicitly rejected without moving geometry. The matched replay starts from the exact saved P5 approaching sample at 3.2328727626623765 nm. It preserves the old 0.02354 nm historical trajectory as a separate comparison. The 20.032 ms validation window includes renewed approach; 1.5/2/3 nm use the same window and dt/dt2/dt4. Handoff events depend on timestep; event/time-accuracy convergence is NOT_ESTABLISHED. All scientific sensitivity and manual visual acceptance remain for user review. No production timestep, neighbor cutoff/skin, glycocalyx or nonspherical lubrication is frozen.

See [中文审核报告](reports/particle6_5/PARTICLE6_5_REVIEW.md) and [machine evidence](reports/particle6_5/PARTICLE6_5_VALIDATION.json).
