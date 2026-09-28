"""Map an accepted CFD state to a new mesh ONLY as an explicitly recorded initial guess."""
from pathlib import Path
import argparse,json,hashlib,shutil,time,xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,d):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def interpolate(source,mesh_path):
 source=Path(source);a=np.load(source/'SV_MESH/mesh_arrays.npz');f=np.load(source/'frozen_flow/flow_arrays_si.npz');target=np.load(mesh_path);x=a['points_m'];t=a['tetra'];y=target['points_m'];u=f['velocity_m_s'];p=f['pressure_pa'];assert np.array_equal(x,f['points_m']) and np.array_equal(t,f['tetra'])
 grid=pv.UnstructuredGrid(np.column_stack([np.full(len(t),4),t]).ravel(),np.full(len(t),pv.CellType.TETRA,np.uint8),x)
 owner=grid.find_containing_cell(y);outside=owner<0;evaluation=y.copy();dist=np.zeros(len(y))
 if outside.any():
  ids,closest=grid.find_closest_cell(y[outside],return_closest_point=True);owner[outside]=ids;evaluation[outside]=closest;dist[outside]=np.linalg.norm(y[outside]-closest,axis=1)
 tet=t[owner];xx=x[tet];A=(xx[:,1:]-xx[:,:1]).swapaxes(1,2);beta=np.linalg.solve(A,(evaluation-xx[:,0])[...,None])[...,0];weights=np.column_stack([1-beta.sum(axis=1),beta]);assert weights.min()>-1e-7 and weights.max()<1+1e-7
 ug=np.einsum('ni,nij->nj',weights,u[tet]);pg=np.einsum('ni,ni->n',weights,p[tet]);assert np.isfinite(ug).all() and np.isfinite(pg).all()
 return ug,pg,dict(nodes=len(y),inside_source_tetra=int((~outside).sum()),projected_to_closest_source_cell=int(outside.sum()),projection_distance_max_um=float(dist.max()*1e6),projection_distance_p95_um=float(np.quantile(dist,.95)*1e6),barycentric_min=float(weights.min()),barycentric_max=float(weights.max()),method='P1 nodal interpolation inside source tets; outside geometric sliver projected to closest source cell; INITIAL GUESS ONLY, not refined CFD')
def prepare(source,case):
 if (Path(case).resolve()/'reports/user_cancellation.json').exists():raise RuntimeError('SKIPPED_BY_USER: no initialization of cancelled case')
 source=Path(source).resolve();case=Path(case).resolve();assert not (case/'run/solver.log').exists();assert json.loads((source/'reports/execution.json').read_text())['status']=='PASS'
 quality=json.loads((source/'reports/flow_quality.json').read_text());assert quality['accepted_final_and_log_checks']
 manifest=json.loads((source/'frozen_flow/manifest.json').read_text())
 for name,item in manifest['files'].items():assert sha(source/'frozen_flow'/name)==item['sha256']
 for c in [source,case]:
  for rel,h in json.loads((c/'input_hashes.json').read_text()).items():assert sha(c/rel)==h
 a=ET.parse(source/'run/solver.xml');b=ET.parse(case/'run/solver.xml')
 for expr in ['.//Time_step_size','.//Density','.//Viscosity/Value','.//Add_BC[@name="INLET"]/Value','.//Add_BC[@name="OUTLET_01"]/Value','.//Add_BC[@name="OUTLET_02"]/Value','.//Add_BC[@name="OUTLET_03"]/Value']:
  assert a.find(expr).text.strip()==b.find(expr).text.strip(),expr
 archive=case/'reports/zero_initial_state_inputs';assert not archive.exists();archive.mkdir()
 for rel in ['policy.json','input_hashes.json','run/solver.xml']:
  dst=archive/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(case/rel,dst)
 ug,pg,info=interpolate(source,case/'SV_MESH/mesh_arrays.npz');m=np.load(case/'SV_MESH/mesh_arrays.npz');wall=np.unique(m['boundary_triangles'][m['facet_tags']==1]);before=ug[wall].copy();ug[wall]=0.
 g=pv.read(case/'SV_MESH/mesh-complete.mesh.vtu');assert np.array_equal(g.points,m['points_m']);assert np.array_equal(g.cells_dict[pv.CellType.TETRA],m['tetra']);g.point_data['Velocity']=ug;g.point_data['Pressure']=pg;dest=case/'initial_state/flow_guess.vtu';dest.parent.mkdir();g.save(dest)
 mesh=b.getroot().find('Add_mesh')
 for tag in ['Initial_velocities_file_path','Initial_pressures_file_path']:
  assert mesh.find(tag) is None;ET.SubElement(mesh,tag).text='../initial_state/flow_guess.vtu'
 b.write(case/'run/solver.xml',encoding='utf-8',xml_declaration=True)
 # Remove only added initialization tags and verify remaining semantic XML exactly.
 for tag in ['Initial_velocities_file_path','Initial_pressures_file_path']:mesh.remove(mesh.find(tag))
 def signature(node):return (node.tag,tuple(sorted(node.attrib.items())),(node.text or '').strip(),tuple(signature(n) for n in node))
 assert signature(b.getroot())==signature(ET.parse(archive/'run/solver.xml').getroot())
 pol=json.loads((case/'policy.json').read_text());pol['initial_state']='Accepted medium-mesh CFD interpolated as initial guess; wall nodes reset to prescribed no-slip. Full fine CFD still required.';pol['initial_guess_source_case']=source.name;put(case/'policy.json',pol)
 hashes=json.loads((case/'input_hashes.json').read_text());hashes.update({'policy.json':sha(case/'policy.json'),'run/solver.xml':sha(case/'run/solver.xml'),'initial_state/flow_guess.vtu':sha(dest)});put(case/'input_hashes.json',hashes)
 info.update(created_unix=time.time(),source_case=str(source),source_flow_sha256=sha(source/'frozen_flow/flow_arrays_si.npz'),source_XML_sha256=sha(source/'run/solver.xml'),target_mesh_sha256=sha(case/'SV_MESH/mesh_arrays.npz'),initial_guess_sha256=sha(dest),wall_nodes=len(wall),wall_initial_velocity_correction_max_m_s=float(np.linalg.norm(before,axis=1).max()),physical_parameters_BCs_dt_tolerances_unchanged=True,full_fine_CFD_and_independent_postprocessing_still_required=True,units=dict(length='m',velocity='m/s',pressure='Pa'))
 put(case/'reports/initial_guess_provenance.json',info);print(json.dumps(info,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--case',type=Path,required=True);args=p.parse_args();prepare(args.source,args.case)
