import numpy as np
from particle_3d.field import FrozenFEMField
from particle_3d.flowfield_conservation_diagnosis import p1_divergence

def test_particle_p1_divergence_exact(cube,real):
    p,t=cube['points'],cube['tetra'];G=np.array([[3.,2.,-7.],[5.,-4.,2.],[1.,9.,6.]])*100
    u=(p-cube['origin'])@G.T+[.001,.002,-.003]
    f=FrozenFEMField(p,t,u,np.zeros(len(p)));d,_=p1_divergence(p,t,u)
    np.testing.assert_allclose(f.gradients_s_inv,np.broadcast_to(G,f.gradients_s_inv.shape),atol=2e-11,rtol=0)
    np.testing.assert_allclose(d,np.trace(G),atol=2e-11,rtol=0)
    rp,rt=real[-1].points_m,real[-1].tetra
    rd,_=p1_divergence(rp,rt,(rp-rp.mean(axis=0))@G.T)
    np.testing.assert_allclose(rd,np.trace(G),atol=1e-6,rtol=0)
