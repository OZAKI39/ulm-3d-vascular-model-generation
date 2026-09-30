import json
from pathlib import Path
import numpy as np
from pd_clot.geometry import make_cloud
from pd_clot.fragment_fluid import FragmentFluid
from pd_clot.regularization.flow import ConsistentFragmentFluid


def test_compiled_surface_matches_existing_geometry_and_load_mask():
    c=json.loads((Path(__file__).parents[1]/'configs/streaming_regularized_demo_coarse.json').read_text())
    g=make_cloud(c['clot']);rng=np.random.default_rng(442)
    x=g.X+rng.normal(size=g.X.shape)*1e-6;b=rng.uniform(.2,1,len(g.pairs));attached=rng.uniform(size=len(x))>.2
    new=ConsistentFragmentFluid(c);old=FragmentFluid(c);old.traction_provider=new.traction_provider
    a=old.surface_force(g,x,b,g.exposed,attached,.0003);z=new.surface_force(g,x,b,g.exposed,attached,.0003)
    for v,w in zip(a,z):np.testing.assert_allclose(v,w,atol=1e-12,rtol=1e-10)
