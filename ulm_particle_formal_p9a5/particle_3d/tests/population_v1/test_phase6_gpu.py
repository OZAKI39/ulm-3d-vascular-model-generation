import os
from pathlib import Path
import numpy as np
import pytest
torch=pytest.importorskip('torch')
pytestmark=pytest.mark.skipif(not torch.cuda.is_available(),reason='Requires actual CUDA device; never counted as CPU pass')
from particle_3d.particle_shapes import Sphere
from particle_3d.population_v1.pairs import candidates,sphere_gap
from particle_3d.population_v1.gpu_neighbors import neighbors
from particle_3d.population_v1.gpu_state import GPUActivePool
from particle_3d.population_v1.provenance import write_json
from test_phase2 import event,particle,shear_environment


def save(name,result):
    write_json(Path(os.environ.get('POPULATION_VALIDATION_DIR','/tmp/population_v1_validation'))/(name+'.json'),result)


def test_cuda_neighbors_exact_pairs_and_realization_isolation():
    rng=np.random.default_rng(2026092906);n=10000
    xyz=rng.uniform(-50e-6,50e-6,(n,3));r=rng.uniform(.375e-6,2e-6,n);rid=np.arange(n)%4
    cpu=[]
    for realization in range(4):
        ix=np.flatnonzero(rid==realization);shapes={int(i):Sphere(xyz[i],r[i]) for i in ix}
        cpu.extend(candidates(shapes))
    cpu=sorted(cpu)
    x=torch.tensor(xyz,device='cuda');radius=torch.tensor(r,device='cuda');real=torch.tensor(rid,device='cuda')
    pairs,g,diag=neighbors(x,radius,real);got=[tuple(p) for p in pairs.cpu().numpy()]
    assert got==cpu
    expected=np.array([sphere_gap(Sphere(xyz[i],r[i]),Sphere(xyz[j],r[j]),i,j).gap_m for i,j in cpu])
    err=float(np.max(abs(g[:,0].cpu().numpy()-expected),initial=0.))
    assert err<1e-18
    assert np.all(rid[np.array(got)[:,0]]==rid[np.array(got)[:,1]])
    # Batch partitioning by independent realization cannot change within-realization pairs.
    rebuilt=[]
    for realization in range(4):
        ix=torch.nonzero(real==realization).flatten();part,_,_=neighbors(x[ix],radius[ix],real[ix]);mapped=ix[part].cpu().numpy()
        rebuilt.extend(tuple(p) for p in mapped)
    assert sorted(rebuilt)==got
    save('phase6_neighbors',dict(PASS=True,N=n,realizations=4,candidate_count=len(cpu),
        maximum_gap_error_m=err,whole_batch_vs_separate_identical=True,diagnostics=diag,
        scope='SYNTHETIC_GEOMETRY_VALIDATION_NOT_PHYSICAL_PRODUCTION_CONCENTRATION'))


def test_cuda_swept_crossing_and_id_pool():
    x=torch.tensor([[-1e-4,0.,0.],[1e-4,0.,0.]],device='cuda',dtype=torch.float64)
    radius=torch.full((2,),1e-6,device='cuda',dtype=torch.float64);real=torch.zeros(2,device='cuda',dtype=torch.int64)
    delta=torch.tensor([[2e-4,0.,0.],[-2e-4,0.,0.]],device='cuda',dtype=torch.float64)
    pairs,_,_=neighbors(x,radius,real,delta=delta)
    assert pairs.cpu().tolist()==[[0,1]]
    real[1]=1
    assert len(neighbors(x,radius,real,delta=delta)[0])==0
    pool=GPUActivePool(capacity=1);env=shear_environment()
    for k in range(4):
        p=particle(env,k+1,[0.,0.,1e-5]);pool.insert(p,event(p),0.)
    pool.remove_mask(torch.tensor([False,True,False,True],device='cuda'))
    assert pool.particle_id[:pool.n].cpu().tolist()==[1,3]
    assert pool.x.dtype==torch.float64 and pool.x.is_cuda
    with pytest.raises(ValueError,match='STABLE_ID_NEVER_REUSED'):
        p=particle(env,2,[0.,0.,1e-5]);pool.insert(p,event(p),0.)
    save('phase6_state',dict(PASS=True,swept_crossing_found=True,foreign_realization_excluded=True,
        stable_ids_after_removal=[1,3],residence='GPU_DEVICE',precision='float64'))
