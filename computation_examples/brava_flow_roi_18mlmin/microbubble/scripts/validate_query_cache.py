"""Full reference-trajectory and audit-byte check before using exact caches."""
from pathlib import Path
import fcntl,hashlib,json,time
import numpy as np
import campaign
from native_adapter import accelerated_job
from query_cache import install

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    proof=campaign.HERE/'data/query_cache_parity.json'
    identity={name:digest(campaign.HERE/'scripts'/name) for name in ['query_cache.py','native_adapter.py','closest_native.cpp']}
    identity.update(flow_sha256=campaign.FLOW_SHA,cohort_sha256=digest(campaign.HERE/'data/cohort.json'))
    with (campaign.HERE/'data/query_cache.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if proof.exists():
            old=json.loads(proof.read_text())
            if old['identity']==identity and old['PASS']:return
        root,env,protected=campaign.setup()
        data=np.load(campaign.HERE/'pilot/reference/mb_000001/trajectory.npz')['samples']
        ids=np.linspace(0,len(data)-1,256,dtype=int)
        points=data[ids,1:4]
        old_samples=[env.field.sample(p) for p in points]
        old_candidates=[env.wall.candidates(p,1e-5) for p in points]
        info=install(env)
        for p,a,c in zip(points,old_samples,old_candidates):
            b=env.field.sample(p)
            for name,value in vars(a).items():
                other=getattr(b,name)
                if isinstance(value,np.ndarray):assert value.tobytes()==other.tobytes()
                else:assert value==other
            assert c.tobytes()==env.wall.candidates(p,1e-5).tobytes()
            assert env.field.sample(p) is b
        event=json.loads((campaign.HERE/'data/cohort.json').read_text())['events'][0]
        job=accelerated_job(campaign.bind_job(env,dict(role='EXACT_QUERY_CACHE_PARITY',dt_s=.0005,cache_identity=identity)))
        start=time.time()
        row=job((event,str(campaign.HERE/'pilot/query_cache/mb_000001')))
        reference=json.loads((campaign.HERE/'pilot/reference/mb_000001/metrics.json').read_text())
        result=dict(identity=identity,query_positions_checked=len(points),query_values_and_candidate_bytes_equal=True,
            full_trajectory_bytes_equal=row['array_scientific_sha256']==reference['array_scientific_sha256'],
            full_audit_bytes_equal=row['audit_scientific_sha256']==reference['audit_scientific_sha256'],
            wall_s=time.time()-start,cache_info=info()._asdict())
        result['PASS']=result['full_trajectory_bytes_equal'] and result['full_audit_bytes_equal']
        campaign.dump(proof,result)
        print(json.dumps(result),flush=True)
        assert result['PASS']

if __name__=='__main__':main()
