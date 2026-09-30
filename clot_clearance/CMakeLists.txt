cmake_minimum_required(VERSION 3.20)
project(clot_clearance LANGUAGES NONE)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(RUN_ENV "PYTHONPATH=${CMAKE_SOURCE_DIR}" "PYTHONDONTWRITEBYTECODE=1"
    "NUMBA_CACHE_DIR=${CMAKE_BINARY_DIR}/numba_cache" "MPLCONFIGDIR=${CMAKE_BINARY_DIR}/matplotlib"
    "OPENBLAS_NUM_THREADS=1" "OMP_NUM_THREADS=1")
add_custom_target(pd_clot ALL
  COMMAND ${CMAKE_COMMAND} -E env ${RUN_ENV} ${Python3_EXECUTABLE} -B -m pd_clot.build
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} COMMENT "Compile and validate local Numba mechanics")
enable_testing()
add_test(NAME prototype_verification COMMAND ${CMAKE_COMMAND} -E env ${RUN_ENV}
  ${Python3_EXECUTABLE} -B -m pytest -q -p no:cacheprovider --basetemp=${CMAKE_BINARY_DIR}/pytest_tmp ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(prototype_verification PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
