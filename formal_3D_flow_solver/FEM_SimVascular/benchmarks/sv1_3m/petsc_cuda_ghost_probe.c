/* Solver-equivalent creation/options chain. Every accessed value is checked.
 * One-rank solver has no remote ghosts; 2 ranks have one block of real ghosts.
 * CPU/GPU options are runtime choices; no backend-conversion workaround. */
#include <petscvec.h>
static PetscErrorCode Probe(void)
{
  Vec g,l=NULL; PetscMPIInt rank,size; PetscInt bs=2,n=4,ng,ghost,rstart,ls=-1;
  PetscScalar *a; const PetscScalar *v; VecType type,lt="NULL"; PetscBool seq,mpi;
  PetscInt errors=0,total=0; PetscReal forward=0,reverse=0;
  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCheck(size==1||size==2,PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Only ranks 1 or 2 supported by this test");
  ng=size>1?1:0; ghost=(1-rank)*2;
  PetscCall(VecCreateGhostBlock(PETSC_COMM_WORLD,bs,n,PETSC_DECIDE,ng,&ghost,&g));
  PetscCall(VecSetFromOptions(g));
  PetscCall(VecGetType(g,&type));
  PetscCall(PetscObjectBaseTypeCompare((PetscObject)g,VECSEQ,&seq));
  PetscCall(PetscObjectBaseTypeCompare((PetscObject)g,VECMPI,&mpi));
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"TYPE rank=%d type=%s base=%s owned=%d ghost_blocks=%d ghost_index=%d\n",rank,type,seq?"seq":mpi?"mpi":"unknown",(int)n,(int)ng,ng?(int)ghost:-1));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));
  PetscCall(VecGetOwnershipRange(g,&rstart,NULL));
  PetscCall(VecGetArray(g,&a));
  for(PetscInt i=0;i<n;i++)a[i]=rstart+i+1;
  PetscCall(VecRestoreArray(g,&a));
  PetscCall(VecGhostGetLocalForm(g,&l));
  if(l){PetscCall(VecGetLocalSize(l,&ls));PetscCall(VecGetType(l,&lt));}
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"LOCAL rank=%d type=%s size=%d expected=%d present=%d\n",rank,lt,(int)ls,(int)(n+bs*ng),!!l));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));
  /* Exercise the same failing API even when GetLocalForm returned NULL. */
  PetscCall(VecGhostUpdateBegin(g,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(g,INSERT_VALUES,SCATTER_FORWARD));
  PetscCheck(l,PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Missing ghost local form");
  PetscCheck(ls==n+bs*ng,PETSC_COMM_WORLD,PETSC_ERR_ARG_SIZ,"Incorrect owned+ghost local size");
  PetscCall(VecGetArrayRead(l,&v));
  for(PetscInt i=0;i<ls;i++){
    PetscScalar expected=i<n?rstart+i+1:ghost*bs+(i-n)+1;
    PetscReal err=PetscAbsScalar(v[i]-expected);forward=PetscMax(forward,err);if(err>1e-13)errors++;
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"VALUE phase=forward rank=%d local=%d actual=%.17g expected=%.17g\n",rank,(int)i,(double)PetscRealPart(v[i]),(double)PetscRealPart(expected)));
  }
  PetscCall(VecRestoreArrayRead(l,&v));
  PetscCall(VecGetArray(l,&a));
  for(PetscInt c=0;c<bs*ng;c++)a[n+c]=100+10*rank+c;
  PetscCall(VecRestoreArray(l,&a));
  PetscCall(VecGhostRestoreLocalForm(g,&l));
  PetscCall(VecGhostUpdateBegin(g,ADD_VALUES,SCATTER_REVERSE));
  PetscCall(VecGhostUpdateEnd(g,ADD_VALUES,SCATTER_REVERSE));
  PetscCall(VecGhostGetLocalForm(g,&l));
  PetscCall(VecGetArrayRead(l,&v));
  for(PetscInt i=0;i<ls;i++){
    PetscScalar expected=i<n?rstart+i+1:100+10*rank+i-n;
    if(size==2&&i<bs)expected+=100+10*(1-rank)+i;
    PetscReal err=PetscAbsScalar(v[i]-expected);reverse=PetscMax(reverse,err);if(err>1e-13)errors++;
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"VALUE phase=reverse rank=%d local=%d actual=%.17g expected=%.17g\n",rank,(int)i,(double)PetscRealPart(v[i]),(double)PetscRealPart(expected)));
  }
  PetscCall(VecRestoreArrayRead(l,&v));
  PetscCall(VecGhostRestoreLocalForm(g,&l));
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"CHECK rank=%d forward_error=%.17g reverse_error=%.17g local_size=%d expected_size=%d errors=%d\n",rank,(double)forward,(double)reverse,(int)ls,(int)(n+bs*ng),(int)errors));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));
  PetscCallMPI(MPI_Allreduce(&errors,&total,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD));
  PetscCheck(total==0,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Ghost data mismatch");
  PetscCall(VecDestroy(&g));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"GHOST_PROBE PASS\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}
int main(int argc,char **argv){PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));PetscCall(Probe());PetscCall(PetscFinalize());return 0;}
