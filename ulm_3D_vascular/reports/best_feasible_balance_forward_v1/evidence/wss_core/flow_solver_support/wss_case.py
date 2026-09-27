"""Explicit-case SI inputs and provenance for production P1 WSS recovery."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
from . import wss


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def coordinate_identity(serialized, original):
    """Accept only exact original order or the solver's exact float32 writeout."""
    exact = np.array_equal(serialized, original)
    float32 = np.array_equal(serialized, np.asarray(original, dtype=np.float32))
    if not (exact or float32):
        raise ValueError('VTU coordinates/order differ from explicit mesh and exact float32 serialization')
    return dict(exact=exact, exact_float32_serialization=float32,
                max_serialization_difference_m=float(np.max(abs(serialized-original))),
                gradients_use_original_mesh_coordinates=True)


def material(case, config=None):
    case = Path(case).resolve()
    config = Path(config).resolve() if config else case / 'run/solver.xml'
    equation = ET.parse(config).find('.//Add_equation[@type="fluid"]')
    if equation is None:
        raise ValueError('Expected an explicit fluid equation in case XML')
    viscosity = equation.find('Viscosity')
    if viscosity is None or viscosity.get('model') != 'Constant':
        raise ValueError('P1 WSS entry currently requires constant Newtonian viscosity')
    mu = float(viscosity.findtext('Value'))
    rho = float(equation.findtext('Density'))
    if not (np.isfinite(mu) and mu > 0 and np.isfinite(rho) and rho > 0):
        raise ValueError('Invalid SI material constants')
    policy_path = case / 'policy.json'
    policy = json.loads(policy_path.read_text())
    if 'nu_m2_s' in policy and not np.isclose(mu/rho, policy['nu_m2_s'], rtol=1e-12, atol=0):
        raise ValueError('Case XML dynamic viscosity disagrees with case policy nu*rho')
    if 'mu_Pa_s' in policy and not np.isclose(mu, policy['mu_Pa_s'], rtol=1e-12, atol=0):
        raise ValueError('Case XML viscosity disagrees with case policy mu')
    return dict(mu_Pa_s=mu, rho_kg_m3=rho, nu_m2_s=mu/rho,
                config=str(config), config_sha256=sha(config),
                policy=str(policy_path), policy_sha256=sha(policy_path),
                unit_contract='SI: points_m, velocity_m_s, pressure_pa; XML mu in Pa.s, rho in kg/m3')


def recover(points, tetra, triangles, tags, velocity, pressure, mu):
    wall_ids = np.flatnonzero(tags == 1)
    wall = triangles[wall_ids]
    owners = wss.boundary_owners(tetra, wall)
    centers, area, normal = wss.wall_geometry(points, tetra, wall, owners)
    gradient = wss.p1_gradients(points, tetra, velocity)
    traction = wss.tangential_traction(gradient[owners], normal, mu)
    magnitude = np.linalg.norm(traction, axis=1)
    display, weight = wss.nodal_average(wall, magnitude, area, len(points))
    mesh, ids = wss.surface(points, wall)
    mesh.point_data['WSS_display_Pa'] = display[ids]
    mesh.point_data['Pressure_Pa'] = pressure[ids]
    for key, value in {'WSS_raw_Pa': magnitude, 'Tangential_viscous_traction_Pa': traction,
                       'Outward_normal': normal, 'Area_m2': area, 'Parent_tetra_zero_based': owners,
                       'Global_boundary_facet_zero_based': wall_ids}.items():
        mesh.cell_data[key] = value
    assert np.isfinite(magnitude).all() and np.all(weight[ids] > 0)
    return mesh, dict(wall_ids=wall_ids, owners=owners, centers=centers, area=area,
                      normal=normal, gradient=gradient, traction=traction,
                      magnitude=magnitude, display=display, weight=weight, ids=ids)
