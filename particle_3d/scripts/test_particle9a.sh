#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export PYTHONPATH=particle_3d/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
.venv/bin/python -m pytest -q -ra \
 particle_3d/tests/particle9a_2mmps \
 particle_3d/tests/particle5/test_resistance_core.py \
 particle_3d/tests/particle5/test_contact_and_time.py \
 particle_3d/tests/particle6_5 \
 particle_3d/tests/particle7/test_inlet_flux_integral.py \
 particle_3d/tests/particle7/test_flux_sampler_actual_inlet.py \
 particle_3d/tests/particle7/test_mb_finite_size_inlet_admission.py \
 particle_3d/tests/particle7/test_no_birth_below_nearfield_handoff.py \
 particle_3d/tests/particle7/test_mb_cumulative_accumulator.py \
 particle_3d/tests/particle7/test_mb_event_time_independent_of_dt.py \
 particle_3d/tests/particle7/test_injection_restart_sequence.py \
 particle_3d/tests/particle82/test_compute_provenance.py \
 particle_3d/tests/particle82/test_outlet_classification.py \
 particle_3d/tests/particle82/test_resume_does_not_recompute.py \
 particle_3d/tests/particle8_2a/test_unchanged_integration_body.py \
 particle_3d/tests/particle8_2a/test_method_B_preserves_anchor_position.py \
 particle_3d/tests/particle8_2a/test_method_B_size_conditioning.py \
 particle_3d/tests/particle8_2a/test_common_ledger_same_across_methods.py \
 -k 'not plot_generation and not review_artifacts and not handoff_scale_sensitivity_artifact and not bridge_parity' \
 --junitxml=particle_3d/reports/particle9a_2mmps/logs/local_tests.xml
