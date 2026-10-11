#!/bin/bash
set -e
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
exec python /workspace/y_clot/project/verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z/run_campaign.py
