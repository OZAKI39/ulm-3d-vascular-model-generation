/* Small device sparse solve and actual options discovery; no vascular physics. */
#include <petscksp.h>
int main(int argc,char **argv) {
 Mat A; Vec x,b,u; KSP ksp; PetscInt n=1024; KSPConvergedReason reason; PetscReal err;
 PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
 PetscCall(MatCreate(PETSC_COMM_WORLD,&A)); PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
 PetscCall(MatSetFromOptions(A)); PetscCall(MatSeqAIJSetPreallocation(A,5,NULL));
 for(PetscInt i=0;i<n;i++) {
  PetscInt cols[5],m=0; PetscScalar v[5];
  cols[m]=i;v[m++]=4.;
  if(i%32>0){cols[m]=i-1;v[m++]=-1.;} if(i%32<31){cols[m]=i+1;v[m++]=-1.;}
  if(i>=32){cols[m]=i-32;v[m++]=-1.;} if(i<n-32){cols[m]=i+32;v[m++]=-1.;}
  PetscCall(MatSetValues(A,1,&i,m,cols,v,INSERT_VALUES));
 }
 PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
 PetscCall(MatCreateVecs(A,&u,&b));PetscCall(VecDuplicate(u,&x));PetscCall(VecSet(u,1.));PetscCall(MatMult(A,u,b));
 PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));PetscCall(KSPSetOperators(ksp,A,A));PetscCall(KSPSetFromOptions(ksp));
 PetscCall(KSPSolve(ksp,b,x));PetscCall(KSPGetConvergedReason(ksp,&reason));
 PetscCall(VecAXPY(x,-1.,u));PetscCall(VecNorm(x,NORM_INFINITY,&err));
 PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO));PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD));PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
 PetscCall(PetscPrintf(PETSC_COMM_WORLD,"SV13P_GPU_SPARSE reason=%d error_inf=%g\n",(int)reason,(double)err));
 PetscCheck(reason>0 && err<1e-7,PETSC_COMM_WORLD,PETSC_ERR_NOT_CONVERGED,"Standalone did not converge");
 PetscCall(KSPDestroy(&ksp));PetscCall(VecDestroy(&x));PetscCall(VecDestroy(&b));PetscCall(VecDestroy(&u));PetscCall(MatDestroy(&A));
 PetscCall(PetscFinalize());return 0;
}
