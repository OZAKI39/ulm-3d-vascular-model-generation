# 本轮实际生产与验证代码摘录

路径相对 FEM_SimVascular；行号和 SHA256 在生成时核对。完整运行命令见 COMMANDS.md。

## solver_support/src/flow_solver_support/wss.py

SHA256 `21168431303441ac2fe96f020bf56f3a6a9aa4be534997e2c2c073c9304c2718`

```python
   1: """P1 wall traction recovery in SI units; no solver/log/case imports.
   2: 
   3: Extracted unchanged algebra from production diagnostics SHA256
   4: efa12e218bcaf99ba07db7194d3b8da79251da709c0a54a97ad6f6179be6e438.
   5: Surface I/O imports PyVista only when called.
   6: """
   7: import numpy as np
   8: 
   9: 
  10: def p1_gradients(points, tetra, velocity):
  11:     x = points[tetra]
  12:     u = velocity[tetra]
  13:     return np.swapaxes(np.linalg.solve(x[:, 1:] - x[:, :1], u[:, 1:] - u[:, :1]), 1, 2)
  14: 
  15: 
  16: def face_keys(faces):
  17:     a = np.ascontiguousarray(np.sort(faces, axis=1), dtype=np.int64)
  18:     return a.view(np.dtype([('a', '<i8'), ('b', '<i8'), ('c', '<i8')])).ravel()
  19: 
  20: 
  21: def boundary_owners(tetra, triangles):
  22:     faces = tetra[:, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]].reshape(-1, 3)
  23:     keys, first, counts = np.unique(face_keys(faces), return_index=True, return_counts=True)
  24:     query = face_keys(triangles)
  25:     index = np.searchsorted(keys, query)
  26:     assert np.all(index < len(keys)), 'Boundary face absent from tetrahedra'
  27:     assert np.array_equal(keys[index], query)
  28:     assert np.all(counts[index] == 1), 'Boundary triangles must have exactly one adjacent element'
  29:     return first[index] // 4
  30: 
  31: 
  32: def wall_geometry(points, tetra, triangles, owners):
  33:     xyz = points[triangles]
  34:     cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
  35:     norm = np.linalg.norm(cross, axis=1)
  36:     assert np.all(norm > 0)
  37:     normal = cross / norm[:, None]
  38:     center = xyz.mean(axis=1)
  39:     direction = center - points[tetra[owners]].mean(axis=1)
  40:     normal[np.einsum('ij,ij->i', normal, direction) < 0] *= -1
  41:     assert np.all(np.einsum('ij,ij->i', normal, direction) > 0)
  42:     return center, norm / 2, normal
  43: 
  44: 
  45: def tangential_traction(gradient, normal, viscosity):
  46:     """Tangential part of viscous stress times outward fluid normal [Pa]."""
  47:     traction = viscosity * np.einsum('nij,nj->ni', gradient + gradient.swapaxes(1, 2), normal)
  48:     return traction - np.einsum('ij,ij->i', traction, normal)[:, None] * normal
  49: 
  50: 
  51: def nodal_average(triangles, values, areas, npoints):
  52:     weight = np.bincount(triangles.ravel(), weights=np.repeat(areas, 3), minlength=npoints)
  53:     total = np.bincount(triangles.ravel(), weights=np.repeat(areas * values, 3), minlength=npoints)
  54:     result = np.zeros(npoints)
  55:     np.divide(total, weight, out=result, where=weight > 0)
  56:     return result, weight
  57: 
  58: 
  59: def surface(points, triangles):
  60:     import pyvista as pv
  61:     ids, inv = np.unique(triangles, return_inverse=True)
  62:     faces = np.column_stack((np.full(len(triangles), 3), inv.reshape(-1, 3)))
  63:     mesh = pv.PolyData(points[ids], faces)
  64:     mesh.point_data['GlobalNodeID_zero_based'] = ids
  65:     return mesh, ids
```

## solver_support/src/flow_solver_support/wss_case.py

SHA256 `85b97a596e1b39fac1b9b75c496d8a6f1d9339d8a3c77b8267d7dc44e79ef47d`

```python
  14: def coordinate_identity(serialized, original):
  15:     """Accept only exact original order or the solver's exact float32 writeout."""
  16:     exact = np.array_equal(serialized, original)
  17:     float32 = np.array_equal(serialized, np.asarray(original, dtype=np.float32))
  18:     if not (exact or float32):
  19:         raise ValueError('VTU coordinates/order differ from explicit mesh and exact float32 serialization')
  20:     return dict(exact=exact, exact_float32_serialization=float32,
  21:                 max_serialization_difference_m=float(np.max(abs(serialized-original))),
  22:                 gradients_use_original_mesh_coordinates=True)
  23: 
  24: 
  25: def material(case, config=None):
  26:     case = Path(case).resolve()
  27:     config = Path(config).resolve() if config else case / 'run/solver.xml'
  28:     equation = ET.parse(config).find('.//Add_equation[@type="fluid"]')
  29:     if equation is None:
  30:         raise ValueError('Expected an explicit fluid equation in case XML')
  31:     viscosity = equation.find('Viscosity')
  32:     if viscosity is None or viscosity.get('model') != 'Constant':
  33:         raise ValueError('P1 WSS entry currently requires constant Newtonian viscosity')
  34:     mu = float(viscosity.findtext('Value'))
  35:     rho = float(equation.findtext('Density'))
  36:     if not (np.isfinite(mu) and mu > 0 and np.isfinite(rho) and rho > 0):
  37:         raise ValueError('Invalid SI material constants')
  38:     policy_path = case / 'policy.json'
  39:     policy = json.loads(policy_path.read_text())
  40:     if 'nu_m2_s' in policy and not np.isclose(mu/rho, policy['nu_m2_s'], rtol=1e-12, atol=0):
  41:         raise ValueError('Case XML dynamic viscosity disagrees with case policy nu*rho')
  42:     if 'mu_Pa_s' in policy and not np.isclose(mu, policy['mu_Pa_s'], rtol=1e-12, atol=0):
  43:         raise ValueError('Case XML viscosity disagrees with case policy mu')
  44:     return dict(mu_Pa_s=mu, rho_kg_m3=rho, nu_m2_s=mu/rho,
  45:                 config=str(config), config_sha256=sha(config),
  46:                 policy=str(policy_path), policy_sha256=sha(policy_path),
  47:                 unit_contract='SI: points_m, velocity_m_s, pressure_pa; XML mu in Pa.s, rho in kg/m3')
  48: 
  49: 
  50: def recover(points, tetra, triangles, tags, velocity, pressure, mu):
  51:     wall_ids = np.flatnonzero(tags == 1)
  52:     wall = triangles[wall_ids]
  53:     owners = wss.boundary_owners(tetra, wall)
  54:     centers, area, normal = wss.wall_geometry(points, tetra, wall, owners)
  55:     gradient = wss.p1_gradients(points, tetra, velocity)
  56:     traction = wss.tangential_traction(gradient[owners], normal, mu)
  57:     magnitude = np.linalg.norm(traction, axis=1)
  58:     display, weight = wss.nodal_average(wall, magnitude, area, len(points))
  59:     mesh, ids = wss.surface(points, wall)
  60:     mesh.point_data['WSS_display_Pa'] = display[ids]
  61:     mesh.point_data['Pressure_Pa'] = pressure[ids]
  62:     for key, value in {'WSS_raw_Pa': magnitude, 'Tangential_viscous_traction_Pa': traction,
  63:                        'Outward_normal': normal, 'Area_m2': area, 'Parent_tetra_zero_based': owners,
  64:                        'Global_boundary_facet_zero_based': wall_ids}.items():
  65:         mesh.cell_data[key] = value
  66:     assert np.isfinite(magnitude).all() and np.all(weight[ids] > 0)
  67:     return mesh, dict(wall_ids=wall_ids, owners=owners, centers=centers, area=area,
  68:                       normal=normal, gradient=gradient, traction=traction,
  69:                       magnitude=magnitude, display=display, weight=weight, ids=ids)
```

## rotate_visualization/prepare_surface_data.py

SHA256 `abb39b8ce2df7cad1973584e35543356071c7e8d0b99563f8c48f7640d1113ea`

```python
   1: """Recover production P1 WSS from an explicitly identified SI flow and case XML."""
   2: from pathlib import Path
   3: import argparse
   4: import json
   5: import sys
   6: import time
   7: import numpy as np
   8: import pyvista as pv
   9: 
  10: HERE = Path(__file__).resolve().parent
  11: sys.path.insert(0, str(HERE.parent / 'solver_support/src'))
  12: from flow_solver_support import wss
  13: from flow_solver_support.wss_case import material, recover, sha, coordinate_identity
  14: 
  15: 
  16: def main():
  17:     started = time.perf_counter()
  18:     parser = argparse.ArgumentParser(description=__doc__)
  19:     parser.add_argument('--case', type=Path, help='Case owning the solver configuration and policy')
  20:     parser.add_argument('--config', type=Path, help='Default: CASE/run/solver.xml')
  21:     parser.add_argument('--flow', type=Path, help='Explicit solver VTU; default: case frozen VTU')
  22:     parser.add_argument('--arrays', type=Path, help='Default: CASE/frozen_flow/flow_arrays_si.npz')
  23:     parser.add_argument('--mesh', type=Path, help='Default: CASE/SV_MESH/mesh_arrays.npz')
  24:     parser.add_argument('--output', type=Path, default=HERE/'input_data/field_diagnostics')
  25:     args = parser.parse_args()
  26:     context = json.loads((HERE/'BUILD_CONTEXT.json').read_text()) if args.case is None else None
  27:     case = (args.case or Path(context['source_case'])).resolve()
  28:     flow = args.flow
  29:     if flow is None:
  30:         candidates = sorted((case/'frozen_flow').glob('*.vtu'))
  31:         if len(candidates) != 1:
  32:             raise ValueError('Specify --flow when frozen_flow does not contain exactly one VTU')
  33:         flow = candidates[0]
  34:     flow = flow.resolve()
  35:     arrays = (args.arrays or case/'frozen_flow/flow_arrays_si.npz').resolve()
  36:     mesh_path = (args.mesh or case/'SV_MESH/mesh_arrays.npz').resolve()
  37:     props = material(case, args.config)
  38:     manifest_path = case/'frozen_flow/manifest.json'
  39:     if manifest_path.exists() and flow.parent == manifest_path.parent.resolve():
  40:         manifest = json.loads(manifest_path.read_text())
  41:         if 'input_configuration_sha256' in manifest:
  42:             assert props['config_sha256'] == manifest['input_configuration_sha256'], 'Configuration differs from frozen solve'
  43:         if 'mesh_arrays_sha256' in manifest:
  44:             assert sha(mesh_path) == manifest['mesh_arrays_sha256'], 'Mesh differs from frozen solve'
  45:         if 'policy_sha256' in manifest:
  46:             assert props['policy_sha256'] == manifest['policy_sha256'], 'Policy differs from frozen solve'
  47:         for name, entry in manifest['files'].items():
  48:             assert sha(manifest_path.parent/name) == entry['sha256'], name
  49:     input_lock = case/'input_hashes.json'
  50:     if input_lock.exists():
  51:         locked = json.loads(input_lock.read_text())
  52:         for key, path in [('run/solver.xml', Path(props['config'])),
  53:                           ('policy.json', case/'policy.json'), ('SV_MESH/mesh_arrays.npz', mesh_path)]:
  54:             assert key in locked and sha(path) == locked[key], 'Case input lock mismatch: '+key
  55:     if context is not None:
  56:         assert sha(flow) == context['source_field_sha256']
  57:     frozen = np.load(arrays)
  58:     points, tetra, triangles, tags, velocity, pressure = [frozen[k] for k in
  59:         ['points_m', 'tetra', 'boundary_triangles', 'facet_tags', 'velocity_m_s', 'pressure_pa']]
  60:     mesh_arrays = np.load(mesh_path)
  61:     for key in ['points_m', 'tetra', 'boundary_triangles', 'facet_tags']:
  62:         assert np.array_equal(frozen[key], mesh_arrays[key]), key
  63:     grid = pv.read(flow)
  64:     coordinate_record = coordinate_identity(grid.points, points)
  65:     assert np.array_equal(np.sort(grid.cells_dict[pv.CellType.TETRA], axis=1), np.sort(tetra, axis=1))
  66:     assert np.array_equal(grid['Velocity'], velocity) and np.array_equal(grid['Pressure'].reshape(-1), pressure)
  67:     assert np.isfinite(velocity).all() and np.isfinite(pressure).all()
  68:     xyz = points[triangles]
  69:     oriented_area = .5*np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
  70:     facet_flux = np.einsum('ij,ij->i', oriented_area, velocity[triangles].mean(axis=1))
  71:     flows = {int(tag): float(facet_flux[tags == tag].sum()) for tag in np.unique(tags)}
  72:     target = float(json.loads((case/'policy.json').read_text())['Q_target_m3_s'])
  73:     assert target > 0 and 4 in flows, 'SI target flow and inlet tag 4 are required'
  74:     epsilon_Q = abs(-flows[4]-target)/target
  75:     epsilon_mass = abs(sum(flows.values()))/target
  76:     assert epsilon_Q <= 1e-6 and epsilon_mass <= 1e-6, (epsilon_Q, epsilon_mass)
  77:     mesh, d = recover(points, tetra, triangles, tags, velocity, pressure, props['mu_Pa_s'])
  78:     ids = d['ids']; magnitude = d['magnitude']; area = d['area']
  79:     assert np.linalg.norm(velocity[ids], axis=1).max() <= 1e-12
  80:     tangent_error = float(abs(np.einsum('ij,ij->i', d['traction'], d['normal'])).max())
  81:     assert tangent_error < max(1., magnitude.max()) * 1e-12
  82:     out = args.output.resolve(); (out/'data').mkdir(parents=True, exist_ok=True)
  83:     mesh.save(out/'data/wall_wss_si.vtp')
  84:     full, pressure_ids = wss.surface(points, triangles)
  85:     full.point_data['Pressure_Pa'] = pressure[pressure_ids]; full.cell_data['Boundary_tag'] = tags
  86:     full.save(out/'data/pressure_surface_si.vtp')
  87:     record = dict(all_pass=True, case=case.name, source_case=str(case), source_field=str(flow),
  88:         source_field_sha256=sha(flow), flow_arrays=str(arrays), flow_arrays_sha256=sha(arrays),
  89:         mesh_arrays=str(mesh_path), mesh_arrays_sha256=sha(mesh_path), material=props,
  90:         coordinate_identity=coordinate_record,
  91:         viscosity_Pa_s=props['mu_Pa_s'], core_implementation=str(Path(wss.__file__).resolve()),
  92:         core_sha256=sha(wss.__file__), entry_sha256=sha(__file__),
  93:         case_adapter_sha256=sha(HERE.parent/'solver_support/src/flow_solver_support/wss_case.py'),
  94:         pressure_Pa=dict(min=float(pressure.min()), max=float(pressure.max())),
  95:         wss_Pa=dict(raw_min=float(magnitude.min()), raw_max=float(magnitude.max()),
  96:             area_weighted_mean=float(np.average(magnitude, weights=area)),
  97:             display_min=float(d['display'][ids].min()), display_max=float(d['display'][ids].max())),
  98:         wss_method='Exact P1 tetra velocity gradients; tangential viscous traction on outward planar WALL facets',
  99:         wss_display='Area-weighted mean of incident facet magnitudes at original wall nodes; geometry unchanged',
 100:         wss_vector_convention='(I-n*n^T)*mu*(grad(u)+grad(u)^T)*n',
 101:         wall_facets=len(magnitude), tangency_max_error_Pa=tangent_error,
 102:         max_wall_speed_m_s=float(np.linalg.norm(velocity[ids], axis=1).max()),
 103:         source_values_clipped=False, mesh_convergence_not_assessed=True,
 104:         measurements=dict(outward_flows_by_tag_m3_s=flows, epsilon_Q=epsilon_Q, epsilon_mass=epsilon_mass),
 105:         outputs_sha256={str(p.relative_to(out)):sha(p) for p in (out/'data').glob('*.vtp')},
 106:         elapsed_seconds=time.perf_counter()-started)
 107:     (out/'COMPUTE_VALIDATION.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
 108:     print(json.dumps(record,indent=2),flush=True)
 109: 
 110: 
 111: if __name__ == '__main__':
 112:     main()
```

## mesh_generate/scripts/generate_tetra_mesh.py

SHA256 `07f914f5fb01634cc621e0ed0bb64ac088609019f5bd63ddac299415d5fc6dd6`

```python
   1: # Executed only by the official SimVascular embedded Python (3.5 compatible).
   2: import json
   3: import os
   4: import time
   5: import sv
   6: import vtk
   7: 
   8: ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
   9: with open(os.path.join(ROOT, 'reports/mesh_and_flow/geometry_import.json')) as stream:
  10:     assert json.load(stream)['status'] == 'PASS'
  11: with open(os.path.join(ROOT, 'configs/mesh_policy.json')) as stream:
  12:     policy = json.load(stream)
  13: factor = float(os.environ.get('MESH_EDGE_FACTOR', '1.0'))
  14: assert factor in (1.0, 0.8)
  15: if factor == 0.8:
  16:     primary = os.path.join(ROOT, 'outputs/mesh_and_flow/mesh_generation/primary')
  17:     assert os.path.isdir(primary), 'The primary generation must be attempted first'
  18:     assert not os.path.exists(os.path.join(primary, 'volume.vtu')), 'No fallback after a generated primary volume'
  19: attempt = 'primary' if factor == 1.0 else 'fallback'
  20: output = os.path.join(ROOT, 'outputs/mesh_and_flow/mesh_generation', attempt)
  21: os.makedirs(output, exist_ok=True)
  22: assert not os.path.exists(os.path.join(output, 'volume.vtu'))
  23: mesher = sv.meshing.TetGen()
  24: mesher.load_model(os.path.join(ROOT, 'outputs/mesh_and_flow/model/imported.vtp'))
  25: mesher.set_walls([1])
  26: assert sorted(mesher.get_model_face_ids()) == [1, 2, 3, 4, 5]
  27: options = sv.meshing.TetGenOptions(global_edge_size=policy['h_ref_m'] * factor,
  28:                                   surface_mesh_flag=True, volume_mesh_flag=True)
  29: options.use_mmg = False
  30: start = time.time()
  31: mesher.generate_mesh(options)
  32: mesh = mesher.get_mesh()
  33: assert mesh.GetNumberOfCells() > 0
  34: mesher.write_mesh(os.path.join(output, 'volume.vtu'))
  35: def write_surface(data, filename):
  36:     writer = vtk.vtkXMLPolyDataWriter()
  37:     writer.SetFileName(os.path.join(output, filename))
  38:     writer.SetInputData(data)
  39:     assert writer.Write() == 1
  40: write_surface(mesher.get_surface(), 'surface.vtp')
  41: for face_id in mesher.get_model_face_ids():
  42:     write_surface(mesher.get_face_polydata(face_id), 'face_' + str(face_id) + '.vtp')
  43: record = {'status': 'PASS', 'attempt': attempt, 'edge_factor': factor,
  44:           'global_edge_size_m': policy['h_ref_m'] * factor,
  45:           'surface_meshing': True, 'volume_meshing': True,
  46:           'points': mesh.GetNumberOfPoints(), 'cells': mesh.GetNumberOfCells(),
  47:           'elapsed_s': time.time() - start, 'generator': 'sv.meshing.TetGen'}
  48: with open(os.path.join(output, 'generation.json'), 'w') as stream:
  49:     json.dump(record, stream, indent=2)
  50: print('TETRA_MESH_GENERATION=' + json.dumps(record))
```

## wss_validation_v2/scripts/case_common.py

SHA256 `82b66e8ac0a62bb8826df74eafe36acab0a75c83a2b02dd71f8440031343eba0`

```python
   1: from pathlib import Path
   2: import sys,json,hashlib,xml.etree.ElementTree as ET
   3: import numpy as np
   4: import pyvista as pv
   5: V=Path(__file__).resolve().parents[1];C=V.parent
   6: sys.path[:0]=[str(C/'solver_support/src'),str(C/'mesh_generate/src')]
   7: from flow_solver_support import wss
   8: from vascular_validation.mesh_diagnostics import tetra_sicn
   9: ROLES={1:'WALL',2:'OUTLET_03',3:'OUTLET_01',4:'INLET',5:'OUTLET_02'}
  10: MU=.00345312;RHO=1056.;DT=1.5040600916882215e-7
  11: 
  12: def vascular_policy(name,kind):
  13:  source=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/policy.json')
  14:  old=json.loads(source.read_text())
  15:  return dict(case=name,kind=kind,dt_s=old['dt_s'],Q_target_m3_s=old['Q_target_m3_s'],mu_Pa_s=MU,nu_m2_s=old['nu_m2_s'],rho_kg_m3=RHO,A_in_m2=old['A_in_m2'],Umean_m_s=old['Umean_m_s'],Dh_m=old['Dh_m'],Re=old['Re'],minimum_steps=20,steady_change_limit=1e-7,steady_consecutive_intervals=3,mass_limit=1e-6,maximum_wall_time_s=43200,maximum_total_steps=500,save_interval_steps=5,MPI_ranks=8,initial_state='zero; no restart or interpolated solution',linear_algebra_backend='CPU PETSc aij/standard; ASM overlap2, subdomain ILU2',dt_rule='Fixed original H0 dt on every mesh; not recalculated from new cap area or mesh size',baseline_policy_sha256=sha(source),baseline_time_step_derivation_reference=old)
  16: 
  17: def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
  18: def dump(p,a):
  19:  p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
  20: def boundary(x,t):
  21:  f=t[:,[[0,1,2],[0,1,3],[0,2,3],[1,2,3]]].reshape(-1,3)
  22:  _,first,count=np.unique(wss.face_keys(f),return_index=True,return_counts=True)
  23:  assert count.max()==2
  24:  b=f[first[count==1]].copy();own=first[count==1]//4
  25:  g=x[b];n=np.cross(g[:,1]-g[:,0],g[:,2]-g[:,0]);direction=g.mean(axis=1)-x[t[own]].mean(axis=1)
  26:  flip=np.einsum('ij,ij->i',n,direction)<0;b[flip]=b[flip][:,[0,2,1]]
  27:  return b,own
  28: 
  29: def write_mesh(case,x,t,b,tags,own):
  30:  case=Path(case);out=case/'SV_MESH';(out/'mesh-surfaces').mkdir(parents=True,exist_ok=True)
  31:  det=np.linalg.det(x[t[:,1:]]-x[t[:,:1]]);assert np.all(det>0)
  32:  q=tetra_sicn(x,t);assert q.min()>0
  33:  grid=pv.UnstructuredGrid(np.column_stack([np.full(len(t),4),t]).ravel(),np.full(len(t),pv.CellType.TETRA,np.uint8),x)
  34:  grid.point_data['GlobalNodeID']=np.arange(1,len(x)+1,dtype=np.int32);grid.cell_data['GlobalElementID']=np.arange(1,len(t)+1,dtype=np.int32);grid.cell_data['minSICN']=q
  35:  grid.save(out/'mesh-complete.mesh.vtu')
  36:  for tag in [None,*np.unique(tags)]:
  37:   sel=np.ones(len(b),bool) if tag is None else tags==tag
  38:   face,ids=wss.surface(x,b[sel]);del face.point_data['GlobalNodeID_zero_based'];face.point_data['GlobalNodeID']=(ids+1).astype(np.int32)
  39:   face.cell_data['GlobalElementID']=(own[sel]+1).astype(np.int32);face.cell_data['ModelFaceID']=tags[sel]
  40:   path=out/'mesh-complete.exterior.vtp' if tag is None else out/'mesh-surfaces'/(ROLES[int(tag)]+'.vtp')
  41:   face.save(path)
  42:  np.savez_compressed(out/'mesh_arrays.npz',points_m=x,tetra=t,boundary_triangles=b,facet_tags=tags,adjacent_tetra=own,min_sicn=q)
  43:  return grid
  44: 
  45: def xml_case(case,dt,steps=500,pipe=False):
  46:  case=Path(case);run=case/'run';run.mkdir(parents=True,exist_ok=True)
  47:  tree=ET.parse(C/'wss_audit/inputs/H0/solver.xml');root=tree.getroot();g=root.find('GeneralSimulationParameters')
  48:  for key,val in {'Time_step_size':dt,'Number_of_time_steps':steps,'Increment_in_saving_VTK_files':5,'Increment_in_saving_restart_files':5}.items():g.find(key).text=str(val)
  49:  if pipe:
  50:   mesh=root.find('Add_mesh')
  51:   for e in list(mesh.findall('Add_face')):
  52:    if e.get('name') not in ['WALL','INLET','OUTLET_03']:mesh.remove(e)
  53:   eq=root.find('Add_equation')
  54:   for e in list(eq.findall('Add_BC')):
  55:    if e.get('name') not in ['WALL','INLET','OUTLET_03']:eq.remove(e)
  56:   inlet=eq.find('Add_BC[@name="INLET"]')
  57:   inlet.find('Profile').text='User_defined';inlet.find('Impose_flux').text='false';inlet.find('Value').text='-1.0'
  58:   ET.SubElement(inlet,'Spatial_profile_file_path').text='inlet_profile.txt'
  59:   outlet=eq.find('Add_BC[@name="OUTLET_03"]')
  60:   for e in list(outlet):outlet.remove(e)
  61:   ET.SubElement(outlet,'Type').text='Trac';ET.SubElement(outlet,'Traction_values_file_path').text='outlet_traction.vtp'
  62:  ET.indent(tree,space='  ');tree.write(run/'solver.xml',encoding='utf-8',xml_declaration=True)
  63:  # Same discretization and KSP/PC as current case; omit redundant debug matrix printing only.
  64:  opts='-ksp_type gmres -ksp_pc_side right -ksp_norm_type unpreconditioned -ksp_rtol 1e-10 -ksp_atol 1e-24 -ksp_max_it 2000 -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_monitor_true_residual -ksp_converged_reason -ksp_view -options_view -options_left -use_gpu_aware_mpi 0 -mat_type aijcusparse -vec_type cuda -ksp_gmres_restart 100 -pc_type asm -pc_asm_overlap 2 -sub_ksp_type preonly -sub_pc_type ilu -sub_pc_factor_levels 2 -sv_pc_adaptive_rebuild true -log_view :petsc_profile.txt -log_view_gpu_time'
  65:  (run/'PETSC_OPTIONS.txt').write_text(opts+'\n')
  66:  if not pipe:
  67:   source=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run/solver.xml')
  68:   def flattened(node,path=''):
  69:    key=path+'/'+node.tag+str(sorted(node.attrib.items()))
  70:    result={key:(node.text or '').strip()} if not len(node) else {}
  71:    for child in node:result.update(flattened(child,key))
  72:    return result
  73:   before=flattened(ET.parse(source).getroot());after=flattened(root)
  74:   changes={k:dict(before=before.get(k),after=after.get(k)) for k in set(before)|set(after) if before.get(k)!=after.get(k)}
  75:   allowed=['Time_step_size','Number_of_time_steps','Increment_in_saving_VTK_files','Increment_in_saving_restart_files']
  76:   assert all('GeneralSimulationParameters' in k and any(k.endswith('/'+s+'[]') for s in allowed) for k in changes),changes
  77:   dump(case/'reports/H0_physical_configuration_identity.json',dict(reference_xml=str(source),reference_sha256=sha(source),new_sha256=sha(run/'solver.xml'),only_time_or_output_changes=True,semantic_changes=changes,material_boundary_equation_and_solver_parameters_identical=True))
  78:  return tree
  79: 
  80: def lock_case(case):
  81:  case=Path(case);paths=[p for root in [case/'SV_MESH',case/'run'] for p in root.rglob('*') if p.is_file()]+[case/'policy.json'];dump(case/'input_hashes.json',{str(p.relative_to(case)):sha(p) for p in paths})
```
