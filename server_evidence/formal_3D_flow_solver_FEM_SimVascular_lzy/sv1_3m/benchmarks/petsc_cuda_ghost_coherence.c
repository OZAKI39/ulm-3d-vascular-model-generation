/* Additional regression: actual device arithmetic, local owned writes, duplicate. */
#include <petscvec.h>
static PetscErrorCode Test(void){
 Vec g,l,d;PetscMPIInt rank;PetscInt ghost,start;PetscScalar *a;const PetscScalar *r;PetscReal norm;
 PetscFunctionBeginUser;PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));ghost=(1-rank)*2;
 PetscCall(VecCreateGhostBlock(PETSC_COMM_WORLD,2,4,PETSC_DECIDE,1,&ghost,&g));PetscCall(VecSetFromOptions(g));
 PetscCall(VecGetOwnershipRange(g,&start,NULL));PetscCall(VecGetArray(g,&a));for(int i=0;i<4;i++)a[i]=start+i+1;PetscCall(VecRestoreArray(g,&a));
 PetscCall(VecScale(g,2)); /* CUDA arithmetic makes device owned values newer. */
 PetscCall(VecGhostUpdateBegin(g,INSERT_VALUES,SCATTER_FORWARD));PetscCall(VecGhostUpdateEnd(g,INSERT_VALUES,SCATTER_FORWARD));
 PetscCall(VecGhostGetLocalForm(g,&l));PetscCall(VecGetArray(l,&a));
 for(int i=0;i<6;i++){
  double expected=2*(i<4?start+i+1:ghost*2+i-4+1);
  PetscCheck(PetscAbsScalar(a[i]-expected)<1e-13,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Device to local mismatch");
  if(i<4)a[i]+=7;else a[i]=10*rank+i-4+1;
 }
 PetscCall(VecRestoreArray(l,&a));PetscCall(VecGhostRestoreLocalForm(g,&l));
 PetscCall(VecNorm(g,NORM_1,&norm));PetscCheck(PetscAbsReal(norm-128)<1e-13,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Local owned writes did not reach device");
 PetscCall(VecGhostUpdateBegin(g,ADD_VALUES,SCATTER_REVERSE));PetscCall(VecGhostUpdateEnd(g,ADD_VALUES,SCATTER_REVERSE));
 PetscCall(VecGhostGetLocalForm(g,&l));PetscCall(VecGetArrayRead(l,&r));
 for(int i=0;i<4;i++){
  double expected=2*(start+i+1)+7+(i<2?10*(1-rank)+i+1:0);
  PetscCheck(PetscAbsScalar(r[i]-expected)<1e-13,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Reverse device owner mismatch");
 }
 PetscCall(VecRestoreArrayRead(l,&r));PetscCall(VecGhostRestoreLocalForm(g,&l));
 PetscCall(VecDuplicate(g,&d));PetscCall(VecCopy(g,d));
 PetscCall(VecGhostUpdateBegin(d,INSERT_VALUES,SCATTER_FORWARD));PetscCall(VecGhostUpdateEnd(d,INSERT_VALUES,SCATTER_FORWARD));
 PetscCall(VecGhostGetLocalForm(d,&l));PetscCall(VecGetArrayRead(l,&r));
 for(int i=0;i<6;i++){
  int gi=i<4?start+i:ghost*2+i-4;int owner=gi/4,li=gi%4;
  double expected=2*(gi+1)+7+(li<2?10*(1-owner)+li+1:0);
  PetscCheck(PetscAbsScalar(r[i]-expected)<1e-13,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Duplicate local mismatch");
 }
 PetscCall(VecRestoreArrayRead(l,&r));PetscCall(VecGhostRestoreLocalForm(d,&l));PetscCall(VecDestroy(&d));PetscCall(VecDestroy(&g));
 PetscCall(PetscPrintf(PETSC_COMM_WORLD,"COHERENCE device_arithmetic=PASS owned_writeback=PASS reverse=PASS duplicate=PASS\n"));
 PetscFunctionReturn(PETSC_SUCCESS);
}
int main(int argc,char**argv){PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));PetscCall(Test());PetscCall(PetscFinalize());return 0;}
