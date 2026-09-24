import numpy as np
import pyvista as pv

def test_export_particle_nodal_identity(real):
    repo,case,mesh,flow,boundaries,field=real
    native=pv.read(case/'run/1-procs/result_071.vtu');export=pv.read(case/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    cp=np.fromfile(case/'run/1-procs/stFile_071.bin',dtype='<f8',offset=56,count=4*mesh.n_points).reshape(-1,4)
    assert np.array_equal(cp[:,:3],native['Velocity'])
    assert np.array_equal(export['Velocity'],field.velocity_nodes_m_s)
    assert np.array_equal(export.points,field.points_m)
    assert np.array_equal(np.sort(export.cells.reshape(-1,5)[:,1:],axis=1),np.sort(field.tetra,axis=1))
    # Independently prescribed barycentric positions exercise Particle sampling.
    for tid in [0,1234,mesh.n_cells-1]:
        w=np.array([.11,.23,.29,.37]);pos=w@field.points_m[field.tetra[tid]]
        value=field.sample(pos)
        assert value.tetra_id==tid
        np.testing.assert_allclose(value.velocity_m_s,w@export['Velocity'][field.tetra[tid]],rtol=1e-10,atol=1e-14)
