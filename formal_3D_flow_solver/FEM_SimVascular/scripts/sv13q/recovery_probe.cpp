// Synthetic fault injection is ONLY in this standalone harness, never in CFD.
#include "sv13q_reuse.h"
#include <cstdio>
struct ProbeContext { void* normal; int failures; };
PetscErrorCode convergence(KSP ksp,PetscInt it,PetscReal norm,KSPConvergedReason* reason,void* ptr) {
  auto* c=static_cast<ProbeContext*>(ptr);
  if (it==0 && c->failures>0) { --c->failures; *reason=KSP_DIVERGED_BREAKDOWN;return PETSC_SUCCESS; }
  return KSPConvergedDefault(ksp,it,norm,reason,c->normal);
}
int main(int argc,char**argv) {
  PetscCall(PetscInitialize(&argc,&argv,nullptr,nullptr));
  Mat A;Vec b,exact,error;KSP ksp;PC pc;
  PetscCall(MatCreateSeqAIJCUSPARSE(PETSC_COMM_SELF,24,24,3,nullptr,&A));
  for(PetscInt i=0;i<24;++i) { PetscCall(MatSetValue(A,i,i,4.,INSERT_VALUES));if(i)PetscCall(MatSetValue(A,i,i-1,-1.,INSERT_VALUES));if(i<23)PetscCall(MatSetValue(A,i,i+1,-.5,INSERT_VALUES)); }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(VecCreateSeqCUDA(PETSC_COMM_SELF,24,&b));PetscCall(VecDuplicate(b,&exact));PetscCall(VecDuplicate(b,&error));PetscCall(VecSet(exact,1.));
  PetscCall(KSPCreate(PETSC_COMM_SELF,&ksp));PetscCall(KSPSetType(ksp,KSPGMRES));PetscCall(KSPGMRESSetRestart(ksp,100));PetscCall(KSPSetPCSide(ksp,PC_RIGHT));PetscCall(KSPSetTolerances(ksp,1e-10,1e-24,PETSC_DEFAULT,2000));PetscCall(KSPSetDiagonalScale(ksp,PETSC_TRUE));PetscCall(KSPSetDiagonalScaleFix(ksp,PETSC_TRUE));PetscCall(KSPGetPC(ksp,&pc));PetscCall(PCSetType(pc,PCASM));PetscCall(PCASMSetOverlap(pc,2));
  PetscCall(PetscOptionsSetValue(nullptr,"-sub_ksp_type","preonly"));PetscCall(PetscOptionsSetValue(nullptr,"-sub_pc_type","ilu"));PetscCall(PetscOptionsSetValue(nullptr,"-sub_pc_factor_levels","2"));PetscCall(KSPSetFromOptions(ksp));
  ProbeContext ctx;ctx.failures=0;PetscCall(KSPConvergedDefaultCreate(&ctx.normal));PetscCall(KSPSetConvergenceTest(ksp,convergence,&ctx,nullptr));
  SV13QReuse state;state.adaptive=PETSC_TRUE;
  PetscCall(MatMult(A,exact,b));PetscCall(sv13q_solve(ksp,A,b,state,1,0));
  PetscCall(MatMult(A,exact,b));ctx.failures=1;PetscCall(sv13q_solve(ksp,A,b,state,2,0));
  KSPConvergedReason reason;PetscCall(KSPGetConvergedReason(ksp,&reason));PetscCall(VecWAXPY(error,-1.,exact,b));PetscReal norm;PetscCall(VecNorm(error,NORM_INFINITY,&norm));
  PetscCheck(reason>0 && norm<1e-9 && state.builds==2 && state.attempts==3,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Recovery probe failed");
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"SV13Q_PROBE recovered=1 solution_error_inf=%.17g builds=%d attempts=%d synthetic_failure=1\n",(double)norm,(int)state.builds,(int)state.attempts));
  PetscCall(MatMult(A,exact,b));ctx.failures=2;PetscCall(sv13q_solve(ksp,A,b,state,3,0));PetscCall(KSPGetConvergedReason(ksp,&reason));
  PetscCheck(reason<0 && state.attempts==5,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Fresh failure must stop after exactly one retry");
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"SV13Q_PROBE fresh_retry_failed_retained=1 no_third_attempt=1 synthetic_failure=1\n"));
  PetscCall(KSPConvergedDefaultDestroy(&ctx.normal));PetscCall(KSPDestroy(&ksp));PetscCall(VecDestroy(&error));PetscCall(VecDestroy(&exact));PetscCall(VecDestroy(&b));PetscCall(MatDestroy(&A));PetscCall(PetscFinalize());return 0;
}
