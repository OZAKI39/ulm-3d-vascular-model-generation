import numpy as np
import pyvista as pv
from sv_validation.geometry import triangles_of
from sv_validation.validation import triangle_flux, normalized_inlet

def test_written_prescription_integral(root, report):
    a=report("inlet_normalization"); v=pv.read(root/a["artifact"])
    actual=-triangle_flux(v.points,triangles_of(v),v.point_data["PrescribedVelocity"])
    assert normalized_inlet(actual,a["Q_target_m3_s"]) <= 1e-10

def test_oriented_flux_in_si_units():
    xyz=np.array([[0,0,0],[1e-6,0,0],[0,2e-6,0.]])
    u=np.tile([0,0,3e-4],(3,1))
    assert abs(triangle_flux(xyz,[[0,1,2]],u)/3e-16-1) < 1e-14
    assert abs(triangle_flux(xyz,[[0,2,1]],u)/3e-16+1) < 1e-14
