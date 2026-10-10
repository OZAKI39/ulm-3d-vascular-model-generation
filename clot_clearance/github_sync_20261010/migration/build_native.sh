#!/bin/bash
set -euo pipefail
source /workspace/y_clot/env.sh
exec > >(tee -a "$Y_CLOT_HOME/logs/native_build.log") 2>&1
trap 'rc=$?; /usr/bin/python3 -c "import json,sys,time;json.dump(dict(exit_code=int(sys.argv[1]),epoch=time.time()),open(sys.argv[2],\"w\"))" "$rc" "$Y_CLOT_HOME/validation/NATIVE_BUILD_EXIT.json"' EXIT
cd "$Y_CLOT_HOME/deps/petsc-source"
export PETSC_DIR="$PWD"
/usr/bin/python3 configure \
  --prefix="$Y_CLOT_HOME/deps/petsc" PETSC_ARCH=arch-y-clot-cuda126 \
  --with-cc="$Y_CLOT_HOME/deps/openmpi/bin/mpicc" \
  --with-cxx="$Y_CLOT_HOME/deps/openmpi/bin/mpicxx" \
  --with-mpiexec="$Y_CLOT_HOME/deps/openmpi/bin/mpiexec" \
  --with-fc=0 --with-fortran-bindings=0 --with-debugging=0 \
  --with-shared-libraries=1 --with-scalar-type=real --with-precision=double \
  --with-64-bit-indices=0 --with-blaslapack-lib='[/usr/lib/x86_64-linux-gnu/liblapack.so.3,/usr/lib/x86_64-linux-gnu/libblas.so.3]' \
  --with-x=0 --with-cxx-dialect=C++17 --with-hip=0 --with-opencl=0 \
  --with-cuda=1 --with-cuda-dir=/usr/local/cuda --with-cudac=/usr/local/cuda/bin/nvcc \
  --with-cuda-arch=89 --with-cuda-dialect=C++17 CUDAFLAGS='-ccbin /usr/bin/g++'
make -j8 PETSC_DIR="$PWD" PETSC_ARCH=arch-y-clot-cuda126 all
make PETSC_DIR="$PWD" PETSC_ARCH=arch-y-clot-cuda126 install
make PETSC_DIR="$Y_CLOT_HOME/deps/petsc" PETSC_ARCH= PETSC_TEST_OPTIONS="-use_gpu_aware_mpi 0" check
printf 'PETSC_BUILD_AND_SELFTEST_PASS\n'
cmake -S "$Y_CLOT_HOME/deps/svMultiPhysics-source" -B "$Y_CLOT_HOME/deps/svMultiPhysics-build" \
  -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 \
  -DSV_LOCAL_VTK_PATH="$Y_CLOT_HOME/deps/solver-inputs/vtk" \
  -DSV_PETSC_DIR="$Y_CLOT_HOME/deps/petsc" \
  "-DSV_ADDITIONAL_CMAKE_ARGS=-DSV_USE_TETGEN:BOOL=OFF;-DTETGEN_LIBRARY_NAME:FILEPATH=$Y_CLOT_HOME/deps/solver-inputs/tetgen/lib/lib_simvascular_thirdparty_tetgen.a;-DCMAKE_CXX_FLAGS:STRING=-I$Y_CLOT_HOME/deps/solver-inputs/tetgen/include;-DBLAS_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/libblas.so.3;-DLAPACK_LIBRARIES:FILEPATH=/usr/lib/x86_64-linux-gnu/liblapack.so.3;-DMPI_C_COMPILER:FILEPATH=$Y_CLOT_HOME/deps/openmpi/bin/mpicc;-DMPI_CXX_COMPILER:FILEPATH=$Y_CLOT_HOME/deps/openmpi/bin/mpicxx"
cmake --build "$Y_CLOT_HOME/deps/svMultiPhysics-build" --parallel 8
printf 'CFD_SOLVER_BUILD_PASS\n'
