#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1 LP_NUM_THREADS=1
export PARTICLE82_RESULTS_ROOT=/root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit PARTICLE82_POINT_ROOT=/root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000
cd /root/particle8_2_runs/cd68cd42f7f7c5a9be4d1da501d0389180de6bb0/checkpoint_pipeline/repo
echo WAITING_FOR_ISOLATED_BENCHMARK
while [ ! -f /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/scaling/SERVER_SCALING_BENCHMARK.json ]; do sleep 10; done
echo 'START continuation'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/run_particle82_continuation.py --previous /root/particle8_2_runs/readonly_inputs/p81/particle_3d/outputs/particle8_1 --output-dir /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/continuation --host-provenance /root/particle8_2_runs/cd68cd42f7f7c5a9be4d1da501d0389180de6bb0/checkpoint_pipeline/host_provenance.json --admission-audit /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/admission/ADMISSION_BASIN_AUDIT.json --scaling-benchmark /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/scaling/SERVER_SCALING_BENCHMARK.json > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/continuation.log 2>&1
echo 'DONE continuation'
echo 'START timestep'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/run_particle82_timestep.py --previous /root/particle8_2_runs/readonly_inputs/p81/particle_3d/outputs/particle8_1 --output-dir /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/timestep --host-provenance /root/particle8_2_runs/cd68cd42f7f7c5a9be4d1da501d0389180de6bb0/checkpoint_pipeline/host_provenance.json --admission-audit /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/admission/ADMISSION_BASIN_AUDIT.json --scaling-benchmark /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/scaling/SERVER_SCALING_BENCHMARK.json --extended-audit /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/continuation/EXTENDED_GUARD_AUDIT.json > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/timestep.log 2>&1
echo 'DONE timestep'
echo 'START natural'
timeout --signal=INT 10800 /root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/run_particle82_natural.py --output-dir /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/natural --host-provenance /root/particle8_2_runs/cd68cd42f7f7c5a9be4d1da501d0389180de6bb0/checkpoint_pipeline/host_provenance.json --point-basin-summary /root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000/POINT_TRACER_SUMMARY.json --admission-audit /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/admission/ADMISSION_BASIN_AUDIT.json --stop-audit /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/admission/stop_audit/SAFETY_STOP_AUDIT.json --extended-audit /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/continuation/EXTENDED_GUARD_AUDIT.json --timestep-audit /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/timestep/TIMESTEP_AUDIT.json --scaling-benchmark /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/scaling/SERVER_SCALING_BENCHMARK.json --budget /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/PREDECLARED_COMPUTE_BUDGET.json > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/natural.log 2>&1
echo 'DONE natural'
echo 'START point_refinement_audit'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/audit_particle82_point_refinement.py --baseline /root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000 --refined /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/point_refinement_100000 --output /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/POINT_REFINEMENT_AUDIT.json > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/point_refinement_audit.log 2>&1
echo 'DONE point_refinement_audit'
echo 'START export'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/export_particle82_results.py --root /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit --host-provenance /root/particle8_2_runs/cd68cd42f7f7c5a9be4d1da501d0389180de6bb0/checkpoint_pipeline/host_provenance.json > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/export.log 2>&1
echo 'DONE export'
echo 'START figures'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/plot_particle82.py --root /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit --baseline-points /root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000 > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/figures.log 2>&1
echo 'DONE figures'
echo 'START animations'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/render_particle82.py --root /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/animations.log 2>&1
echo 'DONE animations'
echo 'START media_audit'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/audit_particle82_media.py --root /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/media_audit.log 2>&1
echo 'DONE media_audit'
echo 'START render_determinism'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/audit_particle82_render_determinism.py --root /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/render_determinism.log 2>&1
echo 'DONE render_determinism'
echo 'START permanent_tests'
/root/particle8_2_runs/env/bin/python -m pytest particle_3d/tests/particle82 -q --junitxml=/root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/particle82_tests.xml > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/permanent_tests.log 2>&1
echo 'DONE permanent_tests'
echo 'START raw_inventory'
/root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/inventory_particle82_raw.py --root /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit --points /root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000 > /root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit/raw_inventory.log 2>&1
echo 'DONE raw_inventory'
echo READY_FOR_ACTUAL_VISUAL_INSPECTION_AND_FINAL_REPORT
