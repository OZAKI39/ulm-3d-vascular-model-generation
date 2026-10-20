#!/bin/bash
set -euo pipefail
/venv/main/bin/python -m venv /workspace/y_clot/envs/inverse_opt_v1
echo /venv/main/lib/python3.12/site-packages > /workspace/y_clot/envs/inverse_opt_v1/lib/python3.12/site-packages/readonly_torch.pth
/workspace/y_clot/envs/inverse_opt_v1/bin/python -m pip install optuna optuna-integration botorch gpytorch
/workspace/y_clot/envs/inverse_opt_v1/bin/python -m pip freeze > /workspace/y_clot/project/verification/physics_constrained_load_inverse_design_v1/20261016T144922Z/optimizer_requirements.txt
