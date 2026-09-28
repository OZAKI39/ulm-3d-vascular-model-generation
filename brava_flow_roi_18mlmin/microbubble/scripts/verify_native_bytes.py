"""Strengthen equality check to include signed-zero / byte representations."""
from pathlib import Path
import json,sys
import numpy as np
HERE=Path(__file__).resolve().parents[1]
config=json.loads((HERE/'config.json').read_text())
sys.path.insert(0,str(Path(config['source_root'])/'particle_3d/src'))
from particle_3d.convex_triangle import triangle_closest_many
from native_adapter import native_closest
tri=np.load(HERE/'data/gpu_mesh_input.npz')['wall_triangles']
rng=np.random.default_rng(2026092806);queries=[]
for k in range(300):
    ids=rng.choice(len(tri),64,replace=False);t=tri[ids]
    queries.append((t[0].mean(axis=0)+rng.normal(size=3)*10**rng.uniform(-10,-5),t))
for t in tri[rng.choice(len(tri),40,replace=False)]:
    for p in [*t,t.mean(axis=0),.5*(t[0]+t[1])]:queries.append((p,t[None]))
different=[0,0];numeric_different=[0,0];nonzero_bit_differences=[0,0]
for p,t in queries:
    for k,(a,b) in enumerate(zip(triangle_closest_many(p,t),native_closest(p,t))):
        different[k]+=a.tobytes()!=b.tobytes()
        numeric_different[k]+=not np.array_equal(a,b)
        mismatch=a.view(np.uint64)!=b.view(np.uint64)
        nonzero_bit_differences[k]+=int(np.count_nonzero(mismatch & ((a!=0)|(b!=0))))
a=json.loads((HERE/'pilot/reference/mb_000001/metrics.json').read_text())
b=json.loads((HERE/'pilot/native/mb_000001/metrics.json').read_text())
result=dict(queries=len(queries),byte_different_kernel_arrays=sum(different),
    byte_different_position_arrays=different[0],byte_different_weight_arrays=different[1],
    numerically_different_arrays=numeric_different,nonzero_bit_differences=nonzero_bit_differences,
    full_kernel_bytes_equal=sum(different)==0,
    interpretation='Native acceleration is used only in first_handoff_event, which discards returned weights. Position bytes must agree. Signed zero in unused weights is reported, not hidden.',
    trajectory_array_scientific_hash_equal=a['array_scientific_sha256']==b['array_scientific_sha256'],
    audit_scientific_hash_equal=a['audit_scientific_sha256']==b['audit_scientific_sha256'],
    production_original_parity_record_unchanged=True)
result['PASS_FOR_ACTUAL_HANDOFF_DEPENDENCY']=different[0]==0 and not any(numeric_different) and not any(nonzero_bit_differences) and result['trajectory_array_scientific_hash_equal'] and result['audit_scientific_hash_equal']
(HERE/'data/native_exact_bytes_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
