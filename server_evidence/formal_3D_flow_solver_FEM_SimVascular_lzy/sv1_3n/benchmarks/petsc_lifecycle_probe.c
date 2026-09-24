/* Small lifecycle contract test. Actual svMP lifecycle is separately traced
 * during the single official GPU smoke; this program is not solver evidence. */
#include <petscvec.h>
#include <stdio.h>
static PetscInt calls = 0;
static PetscErrorCode FinalizeOnce(PetscBool *owned)
{
  PetscBool initialized,finalized;
  PetscMPIInt done;
  PetscFunctionBeginUser;
  if (!*owned) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCallMPI(MPI_Finalized(&done));
  PetscCheck(!done,PETSC_COMM_WORLD,PETSC_ERR_ORDER,"MPI already finalized");
  PetscCall(PetscInitialized(&initialized));
  PetscCall(PetscFinalized(&finalized));
  PetscCheck(initialized && !finalized,PETSC_COMM_WORLD,PETSC_ERR_ORDER,"Invalid PETSc lifecycle");
  PetscCall(PetscFinalize());
  calls++;
  *owned=PETSC_FALSE;
  return PETSC_SUCCESS;
}
int main(int argc,char **argv)
{
  Vec v;
  PetscBool owned=PETSC_TRUE;
  PetscCallMPI(MPI_Init(&argc,&argv));
  PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF,4,&v));
  PetscCall(VecSetFromOptions(v));
  PetscCall(VecSet(v,2));
  PetscCall(VecScale(v,3));
  PetscCall(VecDestroy(&v));
  PetscCall(FinalizeOnce(&owned));
  PetscCall(FinalizeOnce(&owned));
  if(calls!=1)return 8;
  PetscCallMPI(MPI_Finalize());
  printf("LIFECYCLE object_destroyed=1 petsc_finalize_calls=%d mpi_finalized_last=1 PASS\n",(int)calls);
  return 0;
}
