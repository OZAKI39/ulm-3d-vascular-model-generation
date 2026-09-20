/* Diagnosis only: compare application-owned MPI with/without PetscFinalize. */
#include <petscvec.h>
int main(int argc,char **argv){
 Vec v;PetscReal norm;PetscBool proper=PETSC_FALSE;
 MPI_Init(&argc,&argv);PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
 PetscCall(PetscOptionsGetBool(NULL,NULL,"-proper_finalize",&proper,NULL));
 PetscCall(VecCreate(PETSC_COMM_WORLD,&v));PetscCall(VecSetSizes(v,32,PETSC_DECIDE));PetscCall(VecSetFromOptions(v));PetscCall(VecSetUp(v));
 PetscCall(VecSet(v,1));PetscCall(VecScale(v,2));PetscCall(VecNorm(v,NORM_2,&norm));PetscCall(VecDestroy(&v));
 PetscCall(PetscPrintf(PETSC_COMM_WORLD,"LIFECYCLE arithmetic_finished norm=%.17g PetscFinalize=%d\n",(double)norm,(int)proper));
 if(proper)PetscCall(PetscFinalize());
 MPI_Finalize();return 0;
}
