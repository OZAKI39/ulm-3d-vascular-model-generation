import difflib,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256
R=ROOT/'reports/sv1_3q';target=ROOT/'external/sv13q/svMultiPhysics-reuse';base=ROOT/'external/sv13p/svMultiPhysics-reuse';name='Code/Source/solver/petsc_impl.cpp'
old=(base/name).read_text();new=old
# Header follows PETSc compilation guard/includes; locate the first actual PETSc function.
needle='void petsc_initialize('
assert new.count(needle)==1
new=new.replace(needle,'#include "sv13q_reuse.h"\n\n'+needle)
needle='    const PetscInt iEq)\n{'
assert new.count(needle)==1
new=new.replace(needle,'    const PetscInt iEq, SV13QReuse* reuse_policy, const PetscInt timestep)\n{')
needle='    KSPSetUp(psol[cEq].ksp);\n    KSPSolve(psol[cEq].ksp, psol[cEq].b, psol[cEq].b);'
assert new.count(needle)==1
new=new.replace(needle,'''    if (reuse_policy && reuse_policy->enabled()) {
      PetscCallAbort(MPI_COMM_WORLD, sv13q_solve(psol[cEq].ksp, psol[cEq].A,
          psol[cEq].b, *reuse_policy, timestep, cEq));
    } else {
      KSPSetUp(psol[cEq].ksp);
      KSPSolve(psol[cEq].ksp, psol[cEq].b, psol[cEq].b);
    }''')
new=new.replace('    PetscBool reuse_within_timestep_ = PETSC_FALSE;','    SV13QReuse reuse_policy_;\n    PetscBool reuse_within_timestep_ = PETSC_FALSE;')
needle='      "-sv_pc_reuse_within_timestep", &reuse_within_timestep_, nullptr));'
new=new.replace(needle,needle+'''
  PetscCallAbort(MPI_COMM_WORLD, PetscOptionsGetInt(nullptr, nullptr,
      "-sv_pc_rebuild_interval", &reuse_policy_.interval, nullptr));
  PetscCallAbort(MPI_COMM_WORLD, PetscOptionsGetBool(nullptr, nullptr,
      "-sv_pc_adaptive_rebuild", &reuse_policy_.adaptive, nullptr));
  if (reuse_policy_.enabled()) {
    if (reuse_within_timestep_ || (reuse_policy_.adaptive && reuse_policy_.interval) ||
        (!reuse_policy_.adaptive && reuse_policy_.interval != 1 && reuse_policy_.interval != 2 &&
         reuse_policy_.interval != 3 && reuse_policy_.interval != 5)) {
      PetscCallAbort(MPI_COMM_WORLD, PETSC_ERR_ARG_WRONG);
    }
  }
''')
new=new.replace('              com_mod.R.data(), lEq.FSILS.RI.mItr, com_mod.dof, com_mod.cEq);','              com_mod.R.data(), lEq.FSILS.RI.mItr, com_mod.dof, com_mod.cEq,\n              &reuse_policy_, com_mod.cTS);')
assert new!=old;(target/name).write_text(new)
header='Code/Source/solver/sv13q_reuse.h';shutil.copyfile(ROOT/'scripts/sv13q/sv13q_reuse.h',target/header)
patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
patch+=''.join(difflib.unified_diff([], (target/header).read_text().splitlines(True),fromfile='/dev/null',tofile='b/'+header))
p=ROOT/'patches/sv1_3q/ilu_rebuild_policy.patch';p.write_text(patch)
before=json.loads((ROOT/'reports/sv1_3p/source_patch.json').read_text())['after'];after={n:sha256(target/n) for n in before};after[header]=sha256(target/header)
assert [n for n in before if before[n]!=after[n]]==[name]
(R/'source_patch.json').write_text(json.dumps(dict(status='PASS',source=str(target.relative_to(ROOT)),base='Frozen Stage P P1',before=before,after=after,changed_files=[name],added_files=[header],patch=str(p.relative_to(ROOT)),patch_sha256=sha256(p),science_unchanged=True,recovery_audit='Retry happens before RCS rescaling, ghost update or export of correction. Original device RHS saved before stale attempt and restored exactly. PETSc diagonal_scale_fix restores matrix on return even on negative reason. KSP initial guess remains zero. Same tolerance, same operator, only fresh preconditioner.'),indent=2)+'\n')
print('Isolated Stage Q policy patch created: existing CPP plus one header.')
