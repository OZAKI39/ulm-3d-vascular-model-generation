"""Read PETSc's documented binary CSR serialization, compare values; no solver."""
from common import *
import numpy as np
# Identity hashes diagnose byte identity. Floating assembly is also independently
# compared with a machine-roundoff bound (128 epsilon), NOT a relaxed solve tolerance.
put('MATRIX_COMPARISON_PROTOCOL.json',dict(purpose='Resolve first-system byte hash discrepancy before any P2 solve',original_hash_gate='FAIL retained separately',scalar_entry_bound='128 * float64 epsilon * max(abs(A_reference))',rhs_entry_bound='RHS must be bitwise equal; additional norm bound128 epsilon',row_scaled_matrix_bound='max_row(max_abs_delta / max_abs_reference_row) <=128 epsilon',CSR_structure='must be bitwise equal',no_scientific_tolerance_change=True))
def load(p):
 with Path(p).open('rb') as f:
  tag,m,n,nz=np.fromfile(f,dtype='>i4',count=4);assert tag==1211216
  row=np.fromfile(f,dtype='>i4',count=m);col=np.fromfile(f,dtype='>i4',count=nz);v=np.fromfile(f,dtype='>f8',count=nz)
  tagb,nb=np.fromfile(f,dtype='>i4',count=2);assert tagb==1211214 and nb==m
  b=np.fromfile(f,dtype='>f8',count=nb);assert len(f.read())==0
 return row,col,v,b
base=load(HERE/'runtime/p1_baseline/run/first_system.petsc');result={}
for name in ['asm0','schur']:
 p=HERE/f'runtime/p1_{name}/run/first_system.petsc'
 if not p.exists():continue
 a=load(p);structure=np.array_equal(base[0],a[0]) and np.array_equal(base[1],a[1]);d={}
 for label,ref,x in zip(['A','b'],base[2:],a[2:]):
  delta=float(np.max(abs(x-ref)));scale=float(np.max(abs(ref)));d[label]=dict(max_abs_difference=delta,max_abs_reference=scale,scaled_difference=delta/scale,changed_entries=int(np.count_nonzero(x!=ref)),roundoff_units=delta/scale/np.finfo(float).eps,PASS=bool(delta<=128*np.finfo(float).eps*scale))
 ptr=np.r_[0,np.cumsum(base[0])];rowmax=np.maximum.reduceat(abs(base[2]),ptr[:-1]);rowdiff=np.maximum.reduceat(abs(a[2]-base[2]),ptr[:-1]);rowrelative=float(np.max(rowdiff/np.maximum(rowmax,np.finfo(float).tiny)));d['A']['max_row_scaled_difference']=rowrelative;d['A']['row_scaled_roundoff_units']=rowrelative/np.finfo(float).eps;d['A']['PASS']=bool(d['A']['PASS'] and rowrelative<=128*np.finfo(float).eps)
 result[name]=dict(structure_bitwise_equal=structure,values=d,PASS=bool(structure and np.array_equal(base[3],a[3]) and all(x['PASS'] for x in d.values())),byte_identical=sha(p)==sha(HERE/'runtime/p1_baseline/run/first_system.petsc'))
put('P1_MATRIX_RHS_COMPARISON.json',result);print(json.dumps(result,indent=2))
