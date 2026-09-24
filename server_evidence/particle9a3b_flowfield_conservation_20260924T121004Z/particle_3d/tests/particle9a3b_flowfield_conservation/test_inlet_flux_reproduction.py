import json
import numpy as np

def test_inlet_flux_reproduction(real):
    repo,case,mesh,flow,boundaries,field=real;s=boundaries['INLET']
    tri=(np.asarray(s['GlobalNodeID'],int)-1)[s.faces.reshape(-1,4)[:,1:]]
    x=field.points_m[tri];av=np.cross(x[:,1]-x[:,0],x[:,2]-x[:,0])/2
    q=abs(np.einsum('ij,ij->',av,field.velocity_nodes_m_s[tri].mean(axis=1)))
    saved=json.loads((repo/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference/new_flow_flux.json').read_text())['integrated_inlet_Q_m3_s']
    assert np.isclose(q,saved,rtol=1e-13,atol=0)
