"""Inspect the unchanged mixed layout without creating a physical near nullspace."""
import json,re,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3p'
s=ROOT/'external/sv13o/svMultiPhysics/Code/Source/solver';p=s/'petsc_impl.cpp';raw=p.read_text();log=(ROOT/'logs/sv1_3p/remote/P1_SMOKE.log').read_text()
assert 'rows=281452, cols=281452, bs=4' in log and 'type: seqaijcusparse' in log
calls=[]
for f in s.rglob('*'):
 if f.is_file() and f.suffix in ('.cpp','.h','.hpp','.c'):
  if re.search(r'MatSetNearNullSpace\s*\(|MatSetNullSpace\s*\(',f.read_text(errors='replace')):calls.append(str(f.relative_to(ROOT)))
assert not calls
assert 'pindx[i] = jj + dof - 1' in raw and 'uindx[ii+j] = jj + j' in raw
out=dict(status='AUDITED',matrix_type='seqaijcusparse',block_size=4,rows=281452,nodes=70363,DOF_layout='Node-interleaved [ux, uy, uz, p]; existing mapping retained.',near_nullspace_presence=False,near_nullspace_evidence='All adopted solver C/C++ sources contain no MatSetNearNullSpace/MatSetNullSpace; no near-nullspace adapter or options supplied.',existing_fieldsplit_code='Read only to identify existing indices; not selected or changed.',physical_near_nullspace_added=False,DOF_reordered=False,source_path=str(p.relative_to(ROOT)),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),suitability='Unproven: monolithic velocity/pressure coupling; evaluate one conservative PCGAMG profile.',profile='PETSc default aggregation; existing block size 4; no custom near nullspace. Default GPU matrix hierarchy from current CUDA matrix.')
(R/'gamg_matrix_audit.json').write_text(json.dumps(out,indent=2)+'\n')
