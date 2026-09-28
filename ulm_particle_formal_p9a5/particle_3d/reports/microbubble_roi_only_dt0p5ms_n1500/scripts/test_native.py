from pathlib import Path
import json,time
import numpy as np
import campaign
from native_adapter import native_closest,accelerated_job

def main():
    root,env,protected=campaign.setup()
    from particle_3d.convex_triangle import triangle_closest_many
    from particle_3d.formal_cohort_p9a5 import digest
    rng=np.random.default_rng(2026092806)
    queries=[]
    tri=env.wall.triangles
    for k in range(300):
        ids=rng.choice(len(tri),64,replace=False);t=tri[ids]
        p=t[0].mean(axis=0)+rng.normal(size=3)*10**rng.uniform(-10,-5)
        queries.append((p,t))
    for t in tri[rng.choice(len(tri),40,replace=False)]:
        for p in [*t,t.mean(axis=0),.5*(t[0]+t[1])]:queries.append((p,t[None]))
    mismatches=0;maxerror=0.
    for p,t in queries:
        a=triangle_closest_many(p,t);b=native_closest(p,t)
        mismatches+=sum(not np.array_equal(x,y) for x,y in zip(a,b))
        maxerror=max(maxerror,*(float(np.max(abs(x-y))) for x,y in zip(a,b)))
    result=dict(query_count=len(queries),mismatched_arrays=mismatches,max_absolute_difference=maxerror,
        kernel_bitwise_pass=mismatches==0,compiler_flags='-O2 -ffp-contract=off -fno-fast-math')
    campaign.dump(campaign.HERE/'data/native_kernel_parity.json',result)
    print('KERNEL',result,flush=True)
    if mismatches:return
    event=json.loads((campaign.HERE/'data/cohort.json').read_text())['events'][0]
    job=accelerated_job(campaign.bind_job(env,dict(role='COMPILED_GEOMETRY_PARITY',dt_s=.0005)))
    start=time.time();row=job((event,str(campaign.HERE/'pilot/native/mb_000001')))
    a=np.load(campaign.HERE/'pilot/reference/mb_000001/trajectory.npz')['samples']
    b=np.load(campaign.HERE/'pilot/native/mb_000001/trajectory.npz')['samples']
    old=json.loads((campaign.HERE/'pilot/reference/mb_000001/metrics.json').read_text())
    result.update(full_trajectory_bitwise_pass=np.array_equal(a,b),samples_reference=len(a),samples_native=len(b),
        audit_scientific_hash_pass=old['audit_scientific_sha256']==row['audit_scientific_sha256'],
        native_wall_s=time.time()-start,reference_profiled_wall_s=old['wall_seconds'],outlet=row['outlet'])
    result['PASS']=result['full_trajectory_bitwise_pass'] and result['audit_scientific_hash_pass']
    campaign.dump(campaign.HERE/'data/native_kernel_parity.json',result)
    print('FULL_TRAJECTORY',result,flush=True)

if __name__=='__main__':main()
