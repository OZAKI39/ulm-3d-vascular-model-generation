#!/bin/bash
/workspace/y_clot/envs/clot/bin/python /workspace/y_clot/project/analysis/final_v3_v1/driver.py
clot_rank_rc=$?
printf '%s\n' "$clot_rank_rc" > /workspace/y_clot/project/verification/lps_final_native_bugfix_full_v3_regression_v1/20261011T211309Z/C1_final/rank_exit_${OMPI_COMM_WORLD_RANK}.txt
exit "$clot_rank_rc"
