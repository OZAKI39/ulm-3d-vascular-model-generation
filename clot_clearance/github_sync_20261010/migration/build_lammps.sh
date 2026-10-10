#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
exec > >(tee -a "$Y_CLOT_HOME/logs/lammps_build.log") 2>&1
trap 'rc=$?; /usr/bin/python3 -c "import json,sys,time;json.dump(dict(exit_code=int(sys.argv[1]),epoch=time.time()),open(sys.argv[2],\"w\"))" "$rc" "$Y_CLOT_HOME/validation/LAMMPS_BUILD_EXIT.json"' EXIT
cmake -S "$Y_CLOT_HOME/project/tools/lammps-contact-onset-v1/cmake" \
 -B "$Y_CLOT_HOME/project/tools/lammps-native-y-wall-shared-v1-build" \
 -DCMAKE_BUILD_TYPE=Release -DBUILD_MPI=ON -DBUILD_SHARED_LIBS=ON \
 -DPKG_BPM=ON -DPKG_MOLECULE=ON -DDOWNLOAD_POTENTIALS=OFF \
 -DMPI_CXX_COMPILER="$Y_CLOT_HOME/deps/openmpi/bin/mpicxx"
cmake --build "$Y_CLOT_HOME/project/tools/lammps-native-y-wall-shared-v1-build" --parallel 8
