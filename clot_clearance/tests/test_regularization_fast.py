import json,time
from pathlib import Path
import numpy as np
from pd_clot.geometry import make_cloud
from pd_clot.fragment_mechanics import prepare_supported_shape,evaluate_supported as original
from pd_clot.regularization.fast_mechanics import evaluate_supported as fast


def test_scalar_copy_matches_original_and_energy_gradient():
    c=json.loads((Path(__file__).parents[1]/'configs/streaming_regularized_demo_coarse.json').read_text());g=make_cloud(c['clot'])
    rng=np.random.default_rng(194);b=rng.uniform(.4,1,len(g.pairs));s=prepare_supported_shape(g,b,c['safety'])
    x=g.X+rng.normal(size=g.X.shape)*1e-6
    a=original(g,x,b,c['material'],c['safety'],s);z=fast(g,x,b,c['material'],c['safety'],s)
    for v,w in zip(a,z):np.testing.assert_allclose(v,w,atol=3e-10,rtol=1e-8)
    direction=rng.normal(size=x.shape);direction/=np.linalg.norm(direction);eps=1e-9
    def E(y):
        m=fast(g,y,b,c['material'],c['safety'],s);return np.dot(g.volume,m[3]+m[4])
    np.testing.assert_allclose((E(x+eps*direction)-E(x-eps*direction))/(2*eps),-np.sum(z[0]*direction),rtol=1e-5,atol=1e-12)
    measurements={}
    for name,func in [('original',original),('scalar_copy',fast)]:
        t=time.perf_counter()
        for _ in range(50):func(g,x,b,c['material'],c['safety'],s)
        measurements[name]=(time.perf_counter()-t)/50
    print('seconds_per_evaluation',measurements)
