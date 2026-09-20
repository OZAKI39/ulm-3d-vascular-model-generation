"""Host local-form coherence, supported by a failing data probe and PETSc mapping API."""
import subprocess,json,difflib,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];repo=ROOT/'external/compat_cuda/petsc-upstream-sv13m';R=ROOT/'reports/sv1_3m';P=ROOT/'patches/sv1_3m'
files=json.loads((ROOT/'configs/sv1_3m/patch_scope_repair_02.json').read_text())['modified_files']
old={f:subprocess.check_output(['git','-C',str(repo),'show','v3.19.6:'+f],text=True) for f in files}
# All scratch data stays inside the sole development root.
with tempfile.TemporaryDirectory(dir=ROOT/'outputs/sv1_3m') as temp:
 temp=Path(temp)
 for f,s in old.items():p=temp/f;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s)
 subprocess.run(['patch','--batch','--fuzz=0','-p1','-i',str(P/'repair_02/petsc319_cuda_ghost_backport.patch')],cwd=temp,check=True,capture_output=True)
 new={f:(temp/f).read_text() for f in files}
f=files[0];anchor='''    *l         = v->localrep;
'''
addition='''    /* A retained CPU localrep shares the host allocation, not CUDA's device
       allocation. Map the current owned data before exposing that local form. */
    if (*l) {
      PetscBool iscuda;
      PetscCall(PetscObjectTypeCompare((PetscObject)g, VECMPICUDA, &iscuda));
      if (iscuda) {
        PetscScalar *array;
        PetscCall(VecGetArray(g, &array));
        PetscCall(VecRestoreArray(g, &array));
      }
    }
'''
assert new[f].count(anchor)==1;new[f]=new[f].replace(anchor,anchor+addition)
anchor='''  if (*l) {
    PetscCall(VecGhostStateSync_Private(g, *l));
    PetscCall(PetscObjectDereference((PetscObject)*l));
'''
addition='''  if (*l) {
    /* Publish writes through the shared CPU local form. GetArrayWrite does
       not copy stale device values over those host writes. */
    if (*l != g) {
      PetscBool iscuda;
      PetscCall(PetscObjectTypeCompare((PetscObject)g, VECMPICUDA, &iscuda));
      if (iscuda) {
        PetscScalar *array;
        PetscCall(VecGetArrayWrite(g, &array));
        PetscCall(VecRestoreArrayWrite(g, &array));
      }
    }
    PetscCall(VecGhostStateSync_Private(g, *l));
    PetscCall(PetscObjectDereference((PetscObject)*l));
'''
assert new[f].count(anchor)==1;new[f]=new[f].replace(anchor,addition)
patch=''.join(''.join(difflib.unified_diff(old[f].splitlines(True),new[f].splitlines(True),fromfile='a/'+f,tofile='b/'+f)) for f in files)
(P/'repair_03').mkdir(exist_ok=True)
for p in (P/'repair_03/petsc319_cuda_ghost_backport.patch',P/'petsc319_cuda_ghost_backport.patch'):p.write_text(patch)
(ROOT/'configs/sv1_3m/patch_scope_repair_03.json').write_text(json.dumps({'modified_files':files},indent=2)+'\n')
rvec=subprocess.check_output(['git','-C',str(repo),'show','v3.19.6:src/vec/vec/interface/rvector.c'],text=True)
chunks=[]
for name in ('VecGetLocalVector','VecRestoreLocalVector','VecGetArray','VecRestoreArray','VecGetArrayWrite','VecRestoreArrayWrite'):
 a=rvec.index('PetscErrorCode '+name+'(');b=rvec.index('\n}',a)+2;chunks.append(rvec[a:b])
(R/'upstream/official_host_mapping_3.19.6.c').write_text('\n\n'.join(chunks))
D={'iteration':'repair_03','trigger':'repair_02 CUDA MPI reverse values: rank0 [1,2] instead of [111,113], rank1 [5,6] instead of [105,107]','root_cause':'CUDA global owner storage and retained CPU localrep host allocation have separate offload validity; state-number sync alone does not copy device owner updates to the shared host view.','evidence':['ghost_probe_repair_02.json','upstream/official_host_mapping_3.19.6.c','https://petsc.org/release/manualpages/Vec/VecGetArray/','https://gitlab.com/petsc/petsc/-/issues/1096'],'upstream_exact_commit':None,'provenance':'Local minimal-reproducer-supported adaptation of official VecGetLocalVector/GetArray/RestoreArrayWrite host mapping semantics; NOT claimed to be a verbatim upstream ghost fix. Authorized by user section 0.','modified_files':files,'modified_functions':['VecGhostGetLocalForm','VecGhostRestoreLocalForm','VecGhostUpdateBegin','VecGhostUpdateEnd','VecSetType','VecConvert_MPI_MPICUDA_inplace'],'added_lines':sum(x.startswith('+') and not x.startswith('+++') for x in patch.splitlines()),'deleted_lines':sum(x.startswith('-') and not x.startswith('---') for x in patch.splitlines()),'scientific_impact':'NONE — compatibility only','before_result':'CUDA MPI forward PASS, reverse FAIL; CPU and CUDA seq PASS','after_result':'PENDING','tests_added':['same elementwise ghost probe, three repeats','CUDA device arithmetic before local form and local owned writeback, with VecDuplicate']}
(R/'repair_03.json').write_text(json.dumps(D,indent=2)+'\n')
with (R/'upstream_ghost_fix_audit.md').open('a') as out:out.write('''
## repair_03 — host local-form coherence (local adaptation)

The second repair removes API errors but **fails actual reverse data**: rank 0 sees [1,2] instead of [111,113]. The source shows a host `seq` localrep sharing MPI's host allocation while CUDA owner updates invalidate it. Ghost state-number synchronization alone does not transfer data. The local reproducer proves this independently of svMultiPhysics.

Official PETSc 3.19.6 `VecGetLocalVector` maps data with `VecGetArray`, and its restore counterpart publishes writes. The third repair adapts these semantics inside ghost Get/Restore for `mpicuda` with a non-null localrep. Get synchronizes owned data to the existing host view; Restore uses write access to publish local writes without overwriting them from device. Global Vec remains `mpicuda`; no Vec type conversion to CPU, ghost-communication bypass, or solver change is used. This preserves the native CPU localrep architecture, with potential transfers explicitly subject to profiling. CPU paths are unchanged.

This coherence addition is **not** represented as a verbatim upstream commit. It is an independently reproduced, local compatibility repair under user section 0, grounded in official source (`upstream/official_host_mapping_3.19.6.c`). The two true upstream commits remain separately attributed. Added device-arithmetic/owned-writeback/duplicate checks prevent accepting merely suppressed errors.
''')
with (P/'PATCH_PROVENANCE.md').open('a') as out:out.write(f'''\n## Cumulative repair_03 patch\n\nrepair_02 is a CUDA MPI-only adaptation of upstream `{ 'ca0374c5dd2e044b693bca7c31c716bd3be32b2a' }` (v3.20.0). It uses the existing 3.19 CUPM initializer to retain host allocation and ghost metadata. repair_03 adds local-reproducer-supported host-view coherence using official GetArray/RestoreArrayWrite semantics; this part is **not an exact upstream cherry-pick**. Full rationale is in `reports/sv1_3m/repair_03.json`.\n\nCumulative scope: {len(files)} files, {len(D['modified_functions'])} functions, +{D['added_lines']}/-{D['deleted_lines']} lines. Modified files: {', '.join(files)}. A header declaration is included. Per-iteration patches are preserved under repair_01/02/03; the top-level patch is cumulative. Modified hunk line locations are recorded in each unified diff. Scientific impact remains **NONE — compatibility only**.\n''')
print(json.dumps(D))
