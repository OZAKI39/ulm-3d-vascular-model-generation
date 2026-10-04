"""Native solver files for a physical-coordinate inlet and full Cauchy outlet."""
import json
import numpy as np
import pyvista as pv
from audit_pipe_fidelity import boundary_profile


def prepare(mesh, run, inputs, mean, mu):
    radius, length = inputs['radius_m'], inputs['length_m']
    inlet = pv.read(mesh / 'mesh-surfaces/INLET.vtp')
    outlet = pv.read(mesh / 'mesh-surfaces/OUTLET_01.vtp')
    entry, velocity = boundary_profile(inlet, radius, mean)
    exit_trace, _ = boundary_profile(outlet, radius, mean)
    assert entry['PASS'] and exit_trace['PASS']
    ids = np.asarray(inlet.point_data['GlobalNodeID'])
    assert len(np.unique(ids)) == inlet.n_points and ids.min() >= 1
    with (run / 'inlet_exact_velocity.dat').open('w') as stream:
        stream.write(f'3 2 {inlet.n_points}\n0\n10000000000\n')
        for node, speed in zip(ids, velocity):
            stream.write(f'{int(node)}\n0 0 {speed:.17g}\n0 0 {speed:.17g}\n')
    xyz = np.asarray(outlet.points, float)
    gradient = np.zeros((len(xyz),3,3))
    gradient[:,2,:2] = -4 * mean * xyz[:,:2] / radius**2
    pressure = (-8 * mu * mean / radius**2) * (xyz[:,2] - length)
    stress = mu * (gradient + gradient.transpose(0,2,1)) - pressure[:,None,None] * np.eye(3)
    traction = stress[:,:,2]  # n=(0,0,+1) on the planar outlet.
    outlet.point_data['Traction'] = traction
    outlet.save(run / 'outlet_exact_traction.vtp')
    tri = outlet.faces.reshape(-1,4)[:,1:]
    bary = np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])
    qx = np.einsum('qi,tij->tqj',bary,xyz[tri])
    interpolated = np.einsum('qi,tij->tqj',bary,traction[tri])
    exact = np.zeros_like(qx)
    exact[...,:2] = -4 * mu * mean * qx[...,:2] / radius**2
    exact[...,2] = 8 * mu * mean * (qx[...,2]-length) / radius**2
    error = float(np.linalg.norm(interpolated-exact)/np.linalg.norm(exact))
    assert error < 1e-12
    receipt = dict(status='PASS', kind='MANUFACTURED_CONTINUUM_POISEUILLE_WITH_SEPARATELY_BOUNDED_POLYGON_GEOMETRY',
                   inlet=entry, outlet_velocity_trace_representation=exit_trace,
                   outlet_velocity_is_prescribed=False, outlet='t=(-p I + mu*(G+G^T))*n, n=+ez, p(z)=dp_dz*(z-L)',
                   traction_interpolation_relative_L2=error,
                   pressure_gauge='Cauchy normal traction fixes the pressure reference p(L)=0',
                   inlet_parser='General: ndof=3, two constant time samples, original GlobalNodeID',
                   impose_flux=False, nodal_flux_scaling_factor=1.,
                   Q_monitor_m3_s=entry['Q_discrete_m3_s'],Q_analytical_m3_s=entry['Q_exact_m3_s'],
                   monitor_reason='Conservation of prescribed discrete trace uses its own integrated Q; analytical Q error remains a separately reported fidelity/calibration metric. No velocity rescaling.',
                   manufactured_equation='steady convection=0, div(u)=0, dp/dz=mu*Laplacian(u_z)=-8*mu*U/R^2')
    (run.parent / 'reports/BOUNDARY_FIDELITY.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def configure_xml(equation):
    import xml.etree.ElementTree as ET
    for bc in equation.findall('Add_BC'):
        role = bc.get('name')
        if role == 'WALL':
            continue
        for child in list(bc):
            bc.remove(child)
        if role == 'INLET':
            fields = dict(Type='Dir', Time_dependence='General',
                          Temporal_and_spatial_values_file_path='inlet_exact_velocity.dat',
                          Profile='Flat', Impose_flux='false', Zero_out_perimeter='true')
        else:
            assert role == 'OUTLET_01'
            fields = dict(Type='Traction', Traction_values_file_path='outlet_exact_traction.vtp', Traction_multiplier='1.0')
        for key,value in fields.items():
            ET.SubElement(bc,key).text = value
