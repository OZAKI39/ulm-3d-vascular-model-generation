/* HYPRE_config.h.  Generated from HYPRE_config.h.in by configure.  */
/* config/HYPRE_config.h.in.  Generated from configure.in by autoheader.  */

#ifndef HYPRE_CONFIG_H
#define HYPRE_CONFIG_H

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
/* #undef HYPRE_FC_FUNC */

/* As HYPRE_FC_FUNC, but for C identifiers containing underscores. */
/* #undef HYPRE_FC_FUNC_ */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if using MLI */
/* #undef HAVE_MLI */

/* Define to 1 if you have the `MPI_Comm_f2c' function. */
#define HAVE_MPI_COMM_F2C 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if using SuperLU */
/* #undef HAVE_SUPERLU */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 for Alpha platforms */
/* #undef HYPRE_ALPHA */

/* Define to 1 if using long long int for HYPRE_Int and HYPRE_BigInt */
/* #undef HYPRE_BIGINT */

/* Feature branch string */
/* #undef HYPRE_BRANCH_NAME */

/* Define to 1 if using complex values */
/* #undef HYPRE_COMPLEX */

/* Define to 1 if in debug mode */
/* #undef HYPRE_DEBUG */

/* Main development branch? */
/* #undef HYPRE_DEVELOP_BRANCH */

/* Develop branch commit number */
/* #undef HYPRE_DEVELOP_NUMBER */

/* Develop branch string */
/* #undef HYPRE_DEVELOP_STRING */

/* Define to 1 if using OpenMP on device [target alloc version] */
/* #undef HYPRE_DEVICE_OPENMP_ALLOC */

/* Define to 1 if strictly checking OpenMP offload directives */
/* #undef HYPRE_DEVICE_OPENMP_CHECK */

/* Define as follows to set the Fortran name mangling scheme: 0 = unspecified;
   1 = no underscores; 2 = one underscore; 3 = two underscores; 4 = caps, no
   underscores; 5 = one underscore before and after */
#define HYPRE_FMANGLE 0

/* BLAS mangling */
#define HYPRE_FMANGLE_BLAS 2

/* LAPACK mangling */
#define HYPRE_FMANGLE_LAPACK 2

/* Define to 1 if an MPI library is found */
#define HYPRE_HAVE_MPI 1

/* Define to 1 if the routine MPI_Comm_f2c is found */
#define HYPRE_HAVE_MPI_COMM_F2C 1

/* Define to 1 for HP platforms */
/* #undef HYPRE_HPPA */

/* Define to 1 for IRIX64 platforms */
/* #undef HYPRE_IRIX64 */

/* Define to 1 for Linux platform */
#define HYPRE_LINUX 1

/* Define to 1 for Linux on platforms running any version of CHAOS */
/* #undef HYPRE_LINUX_CHAOS */

/* Define to 1 if using quad precision values for HYPRE_Real */
/* #undef HYPRE_LONG_DOUBLE */

/* Define to be the max dimension size (must be at least 3) */
#define HYPRE_MAXDIM 3

/* Define to 1 if using long long int for HYPRE_BigInt */
/* #undef HYPRE_MIXEDINT */

/* Define to 1 if using mixed-precision for HYPRE_Real type */
/* #undef HYPRE_MIXED_PRECISION */

/* Print HYPRE errors */
/* #undef HYPRE_PRINT_ERRORS */

/* Bug reports */
#define HYPRE_RELEASE_BUGS "https://github.com/hypre-space/hypre/issues"

/* Date of release */
#define HYPRE_RELEASE_DATE "2026/01/23"

/* Release name */
#define HYPRE_RELEASE_NAME "hypre"

/* Release number */
#define HYPRE_RELEASE_NUMBER 30100

/* Time of release */
#define HYPRE_RELEASE_TIME "00:00:00"

/* Release version */
#define HYPRE_RELEASE_VERSION "3.1.0"

/* Define to 1 for RS6000 platforms */
/* #undef HYPRE_RS6000 */

/* Disable MPI, enable serial codes. */
/* #undef HYPRE_SEQUENTIAL */

/* Define to 1 if using single precision values for HYPRE_Real */
/* #undef HYPRE_SINGLE */

/* Define to 1 for Solaris. */
/* #undef HYPRE_SOLARIS */

/* Define to 1 if want to use host in hypre tests */
/* #undef HYPRE_TEST_USING_HOST */

/* Using HYPRE timing routines */
/* #undef HYPRE_TIMING */

/* Define to 1 if Caliper instrumentation is enabled */
/* #undef HYPRE_USING_CALIPER */

/* Define to 1 if using cuBLAS */
#define HYPRE_USING_CUBLAS 1

/* Define to 1 if executing on device with CUDA */
#define HYPRE_USING_CUDA 1

/* Define to 1 if using streams */
#define HYPRE_USING_CUDA_STREAMS 1

/* Define to 1 if using cuRAND */
#define HYPRE_USING_CURAND 1

/* Define to 1 if using cuSOLVER */
/* #undef HYPRE_USING_CUSOLVER */

/* Define to 1 if using cuSPARSE */
#define HYPRE_USING_CUSPARSE 1

/* Define to 1 if using device async malloc */
/* #undef HYPRE_USING_DEVICE_MALLOC_ASYNC */

/* Define to 1 if using device memory without UM */
#define HYPRE_USING_DEVICE_MEMORY 1

/* Define to 1 if executing on device with OpenMP */
/* #undef HYPRE_USING_DEVICE_OPENMP */

/* Define to 1 if using DSuperLU */
/* #undef HYPRE_USING_DSUPERLU */

/* Using dxml for Blas */
/* #undef HYPRE_USING_DXML */

/* Using ESSL for Lapack */
/* #undef HYPRE_USING_ESSL */

/* Define to 1 if enabling floating-point exception trap */
/* #undef HYPRE_USING_FPE_TRAP */

/* Define to 1 if executing on GPU device */
#define HYPRE_USING_GPU 1

/* Define to 1 if using GPU aware MPI */
/* #undef HYPRE_USING_GPU_AWARE_MPI */

/* HIP being used */
/* #undef HYPRE_USING_HIP */

/* Define to 1 if using hopscotch hashing */
/* #undef HYPRE_USING_HOPSCOTCH */

/* Define to 1 if using host memory only */
/* #undef HYPRE_USING_HOST_MEMORY */

/* Using internal HYPRE routines */
/* #undef HYPRE_USING_HYPRE_BLAS */

/* Using internal HYPRE routines */
/* #undef HYPRE_USING_HYPRE_LAPACK */

/* Define to 1 if executing on host/device with KOKKOS */
/* #undef HYPRE_USING_KOKKOS */

/* Define to 1 if MAGMA is enabled */
/* #undef HYPRE_USING_MAGMA */

/* Define to 1 if want to track memory operations in hypre */
/* #undef HYPRE_USING_MEMORY_TRACKER */

/* Define to 1 if Node Aware MPI library is used */
/* #undef HYPRE_USING_NODE_AWARE_MPI */

/* Define to 1 if using NVIDIA Tools Extension (NVTX) */
/* #undef HYPRE_USING_NVTX */

/* onemkl::BLAS being used */
/* #undef HYPRE_USING_ONEMKLBLAS */

/* onemkl::rng being used */
/* #undef HYPRE_USING_ONEMKLRAND */

/* onemkl::SPARSE being used */
/* #undef HYPRE_USING_ONEMKLSPARSE */

/* Enable OpenMP support */
/* #undef HYPRE_USING_OPENMP */

/* Define to 1 if using persistent communication */
/* #undef HYPRE_USING_PERSISTENT_COMM */

/* Define to 1 if executing on host/device with RAJA */
/* #undef HYPRE_USING_RAJA */

/* rocBLAS being used */
/* #undef HYPRE_USING_ROCBLAS */

/* rocRAND being used */
/* #undef HYPRE_USING_ROCRAND */

/* rocSOLVER being used */
/* #undef HYPRE_USING_ROCSOLVER */

/* rocSPARSE being used */
/* #undef HYPRE_USING_ROCSPARSE */

/* Define to 1 if using AMD rocTX profiling */
/* #undef HYPRE_USING_ROCTX */

/* SYCL being used */
/* #undef HYPRE_USING_SYCL */

/* Define to 1 if using thrust par_nosync policy */
/* #undef HYPRE_USING_THRUST_NOSYNC */

/* Define to 1 if using UMPIRE */
/* #undef HYPRE_USING_UMPIRE */

/* Define to 1 if using UMPIRE for device memory */
/* #undef HYPRE_USING_UMPIRE_DEVICE */

/* Define to 1 if using UMPIRE for host memory */
/* #undef HYPRE_USING_UMPIRE_HOST */

/* Define to 1 if using UMPIRE for pinned memory */
/* #undef HYPRE_USING_UMPIRE_PINNED */

/* Define to 1 if using UMPIRE for unified memory */
/* #undef HYPRE_USING_UMPIRE_UM */

/* Define to 1 if using unified memory */
/* #undef HYPRE_USING_UNIFIED_MEMORY */

/* GPU warp size */
#define HYPRE_WARP_SIZE 32

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "hypre"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "hypre 3.1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "hypre"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "3.1.0"

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

#endif /* HYPRE_CONFIG_H */
