#ifndef SV13Q_REUSE_H
#define SV13Q_REUSE_H
#include <petscksp.h>

// Stage Q changes only the decision to rebuild the existing preconditioner.
struct SV13QReuse {
  PetscInt interval = 0; // 0: disabled, 1/2/3/5: fixed policy
  PetscBool adaptive = PETSC_FALSE;
  PetscInt last_build = -1, reference = 0, logical = 0, builds = 0, attempts = 0;
  PetscBool pending = PETSC_FALSE;
  bool enabled() const { return adaptive || interval > 0; }
  const char* policy() const { return adaptive ? "RA" : interval == 1 ? "R1" : interval == 2 ? "R2" : interval == 3 ? "R3" : "R5"; }
  const char* rebuild_reason(PetscInt step) const {
    if (last_build < 0) return "FIRST_SOLVE";
    if (adaptive && pending) return "ITERATION_GROWTH";
    if (step-last_build >= (adaptive ? 5 : interval)) return adaptive ? "MAX_AGE" : "FIXED_INTERVAL";
    return nullptr;
  }
};

// Keeps original RHS until the stale attempt is accepted. No failed correction is exported.
inline PetscErrorCode sv13q_solve(KSP ksp, Mat matrix, Vec rhs, SV13QReuse& state,
                                PetscInt step, PetscInt equation) {
  PetscFunctionBeginUser;
  const char* why = state.rebuild_reason(step);
  PetscBool reuse = why ? PETSC_FALSE : PETSC_TRUE;
  Vec saved_rhs = nullptr;
  if (state.adaptive && reuse) {
    PetscCall(VecDuplicate(rhs, &saved_rhs));
    PetscCall(VecCopy(rhs, saved_rhs));
  }
  ++state.logical;
  for (PetscInt attempt=0; attempt<2; ++attempt) {
    const PetscInt age = state.last_build < 0 ? 0 : step-state.last_build;
    const PetscInt old_reference = state.reference;
    if (!reuse) { state.last_build=step; ++state.builds; state.pending=PETSC_FALSE; }
    ++state.attempts;
    PC pc=nullptr;
    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(KSPSetReusePreconditioner(ksp,reuse));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "SV13Q_BEGIN logical=%d attempt=%d step=%d equation=%d policy=%s reuse=%d age=%d ref=%d rebuild_reason=%s builds=%d attempts=%d ksp=%p pc=%p\n",
      (int)state.logical,(int)attempt,(int)step,(int)equation,state.policy(),(int)reuse,(int)age,
      (int)old_reference,why ? why : "REUSE",(int)state.builds,(int)state.attempts,(void*)ksp,(void*)pc));
    PetscCall(KSPSetOperators(ksp,matrix,matrix));
    PetscCall(KSPSetUp(ksp));
    PetscCall(KSPSolve(ksp,rhs,rhs));
    PetscInt iterations=0, maximum=0;
    KSPConvergedReason reason;
    PetscCall(KSPGetIterationNumber(ksp,&iterations));
    PetscCall(KSPGetConvergedReason(ksp,&reason));
    PetscCall(KSPGetTolerances(ksp,nullptr,nullptr,nullptr,&maximum));
    const bool healthy = reason > 0 && iterations < maximum;
    if (healthy && !reuse) state.reference=iterations;
    const double ratio = state.reference > 0 ? double(iterations)/state.reference : 0.;
    if (healthy && state.adaptive && reuse && iterations > 1.5*state.reference) state.pending=PETSC_TRUE;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "SV13Q_END logical=%d attempt=%d step=%d reason=%d iterations=%d healthy=%d ref=%d ratio=%.17g pending=%d recovery=%s\n",
      (int)state.logical,(int)attempt,(int)step,(int)reason,(int)iterations,(int)healthy,
      (int)state.reference,ratio,(int)state.pending,
      attempt ? (healthy ? "STALE_ILU_RECOVERED" : "FRESH_RETRY_FAILED") : "NONE"));
    if (healthy || !state.adaptive || !reuse || attempt==1) break;
    PetscCall(VecCopy(saved_rhs,rhs));
    PetscBool equal=PETSC_FALSE;
    PetscCall(VecEqual(saved_rhs,rhs,&equal));
    PetscCheck(equal,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Stage Q original RHS restore failed");
    PetscBool nonzero=PETSC_FALSE;
    PetscCall(KSPGetInitialGuessNonzero(ksp,&nonzero));
    PetscCheck(!nonzero,PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"Stage Q retry requires frozen zero initial guess");
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"SV13Q_RECOVERY logical=%d step=%d old_reason=%d old_iterations=%d original_rhs_restored=1 initial_guess_zero=1\n",
      (int)state.logical,(int)step,(int)reason,(int)iterations));
    why="STALE_FAILURE";reuse=PETSC_FALSE;
  }
  PetscCall(VecDestroy(&saved_rhs));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
