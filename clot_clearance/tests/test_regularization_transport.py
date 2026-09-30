import json
from pathlib import Path
import numpy as np
from pd_clot.geometry import make_cloud
from pd_clot.regularization.flow import FlowSample
from pd_clot.regularization.transport import drag_coefficient,fragment_drag_half_step


class UniformFlow:
    def sample(self,x,t):
        return FlowSample(np.tile([.001,0,0],(len(x),1)),np.zeros(len(x)),np.zeros((len(x),3,3)))


def test_drag_stokes_limit_and_momentum_conservative_distribution():
    c=json.loads((Path(__file__).parents[1]/'configs/streaming_regularized_demo_coarse.json').read_text())
    g=make_cloud(c['clot']);ids=np.flatnonzero(~g.fixed)[:12];cat=np.zeros(len(g.X),int);cat[ids]=1
    v=np.zeros_like(g.X);m=g.volume*c['clot']['density_kg_m3']
    beta,d,Re,regime=drag_coefficient(m[ids].sum(),g.volume[ids].sum(),np.zeros(3),c['pipe']['density_kg_m3'],c['pipe']['viscosity_Pa_s'])
    assert Re==0;np.testing.assert_allclose(beta,3*np.pi*c['pipe']['viscosity_Pa_s']*d)
    result,rows=fragment_drag_half_step(g,g.X,v,m,[dict(particle_ids=ids.tolist(),fragment_id=8,attached_to_base=False)],cat,UniformFlow(),0,1e-5,c)
    assert len(rows)==1 and np.all(result[cat==0]==0)
    impulse=np.sum(m[:,None]*(result-v),axis=0)
    np.testing.assert_allclose(impulse,rows[0]['integrated_drag_impulse_kg_m_s'],atol=1e-25)
    com=np.average(g.X[ids],axis=0,weights=m[ids])
    assert np.linalg.norm(np.cross(g.X[ids]-com,m[ids,None]*result[ids]).sum(axis=0))<1e-24
    assert 0<result[ids,0].min()<.001


def test_moderate_Re_matches_documented_correlation():
    beta,d,Re,_=drag_coefficient(1e-6,1e-9,np.array([.1,0,0]),1000,.001)
    assert .1<Re<1000
    np.testing.assert_allclose(beta/(3*np.pi*.001*d),1+.15*Re**.687)
