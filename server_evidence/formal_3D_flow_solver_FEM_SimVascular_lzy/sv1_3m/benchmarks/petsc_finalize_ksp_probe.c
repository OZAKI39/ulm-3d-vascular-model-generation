/* Diagnostic: external MPI and unclosed KSP resources match solver exit lifecycle. */
/* Small SPD tridiagonal system solved with GMRES; actual GPU types are checked. */
#include <petscksp.h>
#include <string.h>
int main(int argc, char **argv) {
  Mat A; Vec x,b,truth,r; KSP ksp; PC pc;
  PetscInt i,lo,hi,its,n=128; PetscReal residual,bnorm,error;
  KSPConvergedReason reason;
  const char *mt,*vt,*kt,*pt;
  PetscBool proper=PETSC_FALSE;
  MPI_Init(&argc,&argv);
  PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-proper_finalize",&proper,NULL));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetType(A,MATAIJCUSPARSE));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  PetscCall(MatGetOwnershipRange(A,&lo,&hi));
  for(i=lo;i<hi;i++) {
    PetscCall(MatSetValue(A,i,i,4.0,INSERT_VALUES));
    if(i>0) PetscCall(MatSetValue(A,i,i-1,-1.0,INSERT_VALUES));
    if(i+1<n) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,PETSC_DECIDE,n));
  PetscCall(VecSetType(x,VECCUDA));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecDuplicate(x,&b));PetscCall(VecDuplicate(x,&truth));PetscCall(VecDuplicate(x,&r));
  PetscCall(VecSet(truth,1.0));PetscCall(MatMult(A,truth,b));PetscCall(VecSet(x,0.0));
  PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));PetscCall(KSPSetOperators(ksp,A,A));
  PetscCall(KSPSetType(ksp,KSPGMRES));PetscCall(KSPGetPC(ksp,&pc));PetscCall(PCSetType(pc,PCJACOBI));
  PetscCall(KSPSetTolerances(ksp,1e-12,1e-14,PETSC_DEFAULT,1000));PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSolve(ksp,b,x));PetscCall(KSPGetConvergedReason(ksp,&reason));PetscCall(KSPGetIterationNumber(ksp,&its));
  PetscCall(MatMult(A,x,r));PetscCall(VecAXPY(r,-1.0,b));PetscCall(VecNorm(r,NORM_2,&residual));PetscCall(VecNorm(b,NORM_2,&bnorm));
  PetscCall(VecAXPY(x,-1.0,truth));PetscCall(VecNorm(x,NORM_INFINITY,&error));
  PetscCall(MatGetType(A,&mt));PetscCall(VecGetType(x,&vt));PetscCall(KSPGetType(ksp,&kt));PetscCall(PCGetType(pc,&pt));
  PetscCall(KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD));PetscCall(PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"SMOKE Mat=%s Vec=%s KSP=%s PC=%s reason=%d iterations=%" PetscInt_FMT " relative_residual=%.17g error_inf=%.17g\n",mt,vt,kt,pt,(int)reason,its,(double)(residual/bnorm),(double)error));
  PetscCheck(strstr(mt,"cusparse") && strstr(vt,"cuda"),PETSC_COMM_WORLD,PETSC_ERR_SUP,"PETSC_GPU_TYPES_NOT_ACTIVE");
  PetscCheck(reason>0 && residual/bnorm<=1e-10 && error<=1e-10,PETSC_COMM_WORLD,PETSC_ERR_NOT_CONVERGED,"GPU smoke solution failed");
  if(proper) {
  PetscCall(KSPDestroy(&ksp));PetscCall(MatDestroy(&A));PetscCall(VecDestroy(&x));PetscCall(VecDestroy(&b));PetscCall(VecDestroy(&truth));PetscCall(VecDestroy(&r));
  PetscCall(PetscFinalize());
  }
  MPI_Finalize();
  return 0;
}
