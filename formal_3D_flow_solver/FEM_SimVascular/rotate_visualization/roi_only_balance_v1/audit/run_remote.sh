#!/bin/bash
set -euo pipefail
cd /workspace/roi_only_flow_visualization_20260928
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow
/root/particle8_2_runs/env/bin/python -B render_visualization.py > logs/render_visualization.log 2>&1
/root/particle8_2_runs/env/bin/python -B render_surface_fields.py > logs/render_surface_fields.log 2>&1
/root/particle8_2_runs/env/bin/python -B validate_results.py > logs/validate_results.log 2>&1
/root/particle8_2_runs/env/bin/python -B validate_surface_fields.py > logs/validate_surface_fields.log 2>&1
