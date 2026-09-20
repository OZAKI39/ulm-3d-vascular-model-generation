"""Capture official Git history before constructing the minimal backport."""
import subprocess,json,hashlib,difflib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];repo=ROOT/'external/compat_cuda/petsc-upstream-sv13m';R=ROOT/'reports/sv1_3m';E=R/'upstream';E.mkdir(exist_ok=True)
fix='8ba8e283c7d3ba8c2bd73cc674b79590a073ab5a';f='src/vec/vec/impls/mpi/commonmpvec.c'
commands=[['remote','-v'],['show','--format=fuller','--stat',fix],['show',fix,'--',f],['log','v3.25.0','--format=%H %ad %s','--date=iso','-S','PetscObjectBaseTypeCompare','--',f],['blame','-L','70,88',fix,'--',f],['diff','v3.24.6','v3.25.0','--',f],['tag','--contains',fix,'--sort=version:refname']]
records=[]
for i,args in enumerate(commands):
 p=subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True);out=E/f'git_{i:02d}.txt';out.write_text(p.stdout+p.stderr);records.append(dict(command=['git',*args],exit_code=p.returncode,artifact=str(out.relative_to(ROOT)),sha256=hashlib.sha256(out.read_bytes()).hexdigest()));assert p.returncode==0
old=subprocess.check_output(['git','-C',str(repo),'show','v3.19.6:'+f],text=True)
new=subprocess.check_output(['git','-C',str(repo),'show',fix+':'+f],text=True)
(E/'commonmpvec_v3.19.6.c').write_text(old);(E/'commonmpvec_first_fixed.c').write_text(new)
# Backport only the five concrete comparisons changed upstream. No wholesale replacement.
patched=old
for name in ('VecGhostGetLocalForm','VecGhostUpdateBegin','VecGhostUpdateEnd'):
 start=patched.index('PetscErrorCode '+name+'(');end=patched.find('\n}',start)+2
 section=patched[start:end].replace('PetscObjectTypeCompare(', 'PetscObjectBaseTypeCompare(')
 patched=patched[:start]+section+patched[end:]
patch=''.join(difflib.unified_diff(old.splitlines(True),patched.splitlines(True),fromfile='a/'+f,tofile='b/'+f))
P=ROOT/'patches/sv1_3m';(P/'repair_01').mkdir(exist_ok=True)
for p in (P/'petsc319_cuda_ghost_backport.patch',P/'repair_01/petsc319_cuda_ghost_backport.patch'):p.write_text(patch)
ref=json.loads((R/'reference_manifest.json').read_text())
d=dict(status='PASS',official_repository='https://gitlab.com/petsc/petsc.git',base_tag='v3.19.6',base_archive_sha256=ref['PETSc_archive_sha256'],first_fixed_commit=fix,first_fixed_tag='v3.25.0',commands=records,modified_files=[f],modified_functions=['VecGhostGetLocalForm','VecGhostUpdateBegin','VecGhostUpdateEnd'],added_lines=5,deleted_lines=5,repair='repair_01',scope='Exact upstream base-type checks; MPI metadata conversion remains separately tested')
(R/'upstream_audit.json').write_text(json.dumps(d,indent=2)+'\n')
(R/'upstream_ghost_fix_audit.md').write_text(f'''# Official PETSc ghost-vector audit

Official repository: https://gitlab.com/petsc/petsc.git

First base-type dispatch fix: [{fix}](https://gitlab.com/petsc/petsc/-/commit/{fix}), 2026-03-13, first released tag **v3.25.0**. Commit rationale: “Vec: add device support for ghost vectors”; closes official issue [1096](https://gitlab.com/petsc/petsc/-/issues/1096). The issue identifies both strict type matching and destroyed ghost metadata during backend conversion. Full git log, blame, commit diff and tag diff are in `upstream/`, with commands/hashes in `upstream_audit.json`.

Old 3.19.6:
```c
PetscCall(PetscObjectTypeCompare((PetscObject)g, VECMPI, &ismpi));
PetscCall(PetscObjectTypeCompare((PetscObject)g, VECSEQ, &isseq));
```
First fixed source:
```c
PetscCall(PetscObjectBaseTypeCompare((PetscObject)g, VECMPI, &ismpi));
PetscCall(PetscObjectBaseTypeCompare((PetscObject)g, VECSEQ, &isseq));
```

`seqcuda`/`mpicuda` are derived types. Concrete equality misses both; the existing base-aware helper matches their seq/mpi families. The Stage M pre-patch CPU probes pass; CUDA probes record NULL local form and UpdateBegin error.

repair_01 changes only five calls in GetLocalForm, UpdateBegin, UpdateEnd. RestoreLocalForm has no concrete type dispatch; it only synchronizes state and releases a reference, so needs no change. `VecGhostIsLocalForm` is unused in the solver/probe path and is not expanded speculatively. The upstream control-flow style rewrite is omitted because it is semantically identical.

This fixes recognition, but does **not** yet establish that MPI CUDA conversion retains ghost metadata. That is explicitly checked by the two-rank data probe; a second controlled repair requires independent evidence. Sequential solver runs have zero ghost entries; two-rank tests verify nonzero remote ghosts and reverse contributions element by element.
''')
(P/'PATCH_PROVENANCE.md').write_text(f'''# PETSc 3.19.6 + upstream CUDA ghost compatibility backport

Base archive SHA256: `{ref['PETSc_archive_sha256']}`.

repair_01 upstream: `{fix}`, first fixed tag `v3.25.0`.

Reason: derived CUDA Vec types are rejected by concrete ghost dispatch.

Modified file: `{f}`. Functions: GetLocalForm, UpdateBegin, UpdateEnd. Five removed lines / five added lines, exact upstream substitutions. No RestoreLocalForm change.

Scientific impact: **NONE — compatibility only**. svMultiPhysics source and all model inputs remain unchanged. Original PETSc tree/archive are read-only; clean extraction is patched in a separate tree.

Revert: in a disposable clean source tree, `patch -R -p1 < petsc319_cuda_ghost_backport.patch`; then discard all build objects and rebuild. Never apply/revert in historical trees.
''')
print(json.dumps(d,indent=2))
