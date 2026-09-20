"""Narrow CUDA MPI in-place conversion adaptation from official PETSc 3.20."""
import subprocess,json,difflib,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];repo=ROOT/'external/compat_cuda/petsc-upstream-sv13m';R=ROOT/'reports/sv1_3m';P=ROOT/'patches/sv1_3m';fix='ca0374c5dd2e044b693bca7c31c716bd3be32b2a'
def show(tag,f):return subprocess.check_output(['git','-C',str(repo),'show',tag+':'+f],text=True)
files=['src/vec/vec/impls/mpi/commonmpvec.c','include/petsc/private/vecimpl.h','src/vec/vec/interface/vecreg.c','src/vec/vec/impls/mpi/cupm/cuda/vecmpicupm.cu']
old={f:show('v3.19.6',f) for f in files};new=dict(old)
f=files[0]
for name in ('VecGhostGetLocalForm','VecGhostUpdateBegin','VecGhostUpdateEnd'):
 start=new[f].index('PetscErrorCode '+name+'(');end=new[f].find('\n}',start)+2;new[f]=new[f][:start]+new[f][start:end].replace('PetscObjectTypeCompare(','PetscObjectBaseTypeCompare(')+new[f][end:]
f=files[1];anchor='#include <petscvec.h>'
assert anchor in new[f];new[f]=new[f].replace(anchor,anchor+'\n#if defined(PETSC_HAVE_CUDA)\nPETSC_INTERN PetscErrorCode VecConvert_MPI_MPICUDA_inplace(Vec);\n#endif',1)
f=files[2];anchor='  /* same reasons for VECCUDA and VECVIENNACL */'
addition='''  /* CUDA-only subset of upstream ca0374c5: keep MPI ghost storage in place. */
#if defined(PETSC_HAVE_CUDA)
  {
    PetscBool srcMPI, dstCUDA, dstMPICUDA;
    PetscCall(PetscObjectTypeCompare((PetscObject)vec, VECMPI, &srcMPI));
    PetscCall(PetscStrcmp(method, VECCUDA, &dstCUDA));
    PetscCall(PetscStrcmp(method, VECMPICUDA, &dstMPICUDA));
    if (srcMPI && ((size > 1 && dstCUDA) || dstMPICUDA)) {
      PetscCall(VecConvert_MPI_MPICUDA_inplace(vec));
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
#endif
'''
assert anchor in new[f];new[f]=new[f].replace(anchor,addition+anchor,1)
f=files[3];anchor='PetscErrorCode VecCUDAGetArrays_Private('
addition='''/* Upstream ca0374c5 in-place initialization, adapted to the 3.19 CUPM API.
   Retain the existing host allocation, localrep and localupdate. */
PetscErrorCode VecConvert_MPI_MPICUDA_inplace(Vec v)
{
  PetscDeviceContext dctx;
  PetscFunctionBegin;
  PetscCall(PetscDeviceContextGetCurrentContextAssertType_Internal(&dctx, PETSC_DEVICE_CUDA));
  PetscCall(VecMPI_CUDA.Initialize_CUPMBase(v, PETSC_FALSE, static_cast<Vec_MPI *>(v->data)->array, nullptr, dctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

'''
assert anchor in new[f];new[f]=new[f].replace(anchor,addition+anchor,1)
patch=''.join(''.join(difflib.unified_diff(old[f].splitlines(True),new[f].splitlines(True),fromfile='a/'+f,tofile='b/'+f)) for f in files)
(P/'repair_02').mkdir(exist_ok=True);(P/'repair_02/petsc319_cuda_ghost_backport.patch').write_text(patch)
(ROOT/'configs/sv1_3m/patch_scope_repair_02.json').write_text(json.dumps({'modified_files':files},indent=2)+'\n')
commands=[['show','--format=fuller',fix,'--','include/petsc/private/veccupmimpl.h',files[2],files[3]],['tag','--contains',fix,'--sort=version:refname'],['blame','-L','100,115','v3.20.0','--',files[2]]]
for i,args in enumerate(commands):(R/'upstream'/f'conversion_{i}.txt').write_text(subprocess.check_output(['git','-C',str(repo),*args],text=True))
D={'iteration':'repair_02','trigger':'repair_01 CUDA MPI NULL local form / Vector is not ghosted (3/3)','root_cause':'3.19.6 VecSetType destroys Vec_MPI data; localrep/localupdate are lost.','upstream_commit':fix,'first_fixed_tag':'v3.20.0','adaptation':'CUDA MPI conversion only using existing 3.19 Initialize_CUPMBase; retain host allocation in place. No HIP/Kokkos, generic conversion rewrite, pinned reallocations or new solver code.','modified_files':files,'modified_functions':['VecGhostGetLocalForm','VecGhostUpdateBegin','VecGhostUpdateEnd','VecSetType','VecConvert_MPI_MPICUDA_inplace'],'added_lines':sum(x.startswith('+') and not x.startswith('+++') for x in patch.splitlines()),'deleted_lines':sum(x.startswith('-') and not x.startswith('---') for x in patch.splitlines()),'scientific_impact':'NONE — compatibility only','before_result':'CPU seq, CUDA seq, CPU MPI PASS; CUDA MPI FAIL','after_result':'PENDING','tests_added':['same deterministic CPU/CUDA sequential and MPI ghost probes, three repeats']}
(R/'repair_02.json').write_text(json.dumps(D,indent=2)+'\n')
with (R/'upstream_ghost_fix_audit.md').open('a') as out:out.write(f'\n## repair_02 — conversion retains ghost metadata\n\nOfficial [{fix}](https://gitlab.com/petsc/petsc/-/commit/{fix}), first tag v3.20.0, adds in-place CUDA conversion. The complete 388-line upstream change is **not** imported. Only the MPI CUDA route is adapted using the existing 3.19 initializer. It preserves the existing host allocation and Vec_MPI metadata, and sets the CUDA operations/type. Forward/reverse data checks remain mandatory. See `repair_02.json` and `upstream/conversion_*.txt`.\n')
print(json.dumps(D))
