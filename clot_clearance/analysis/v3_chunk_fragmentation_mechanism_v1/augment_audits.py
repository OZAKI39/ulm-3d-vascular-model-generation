from common import *
from mixed_mode import corotated
from graph_ops import topology,boundary
def main():
 r=reference();p=O/'CONTACT_0376MS_AUDIT.json';a=json.loads(p.read_text());lookup={int(s['step']):s for s in [load_npz(p) for p in (O/'derived/recovered_accepted').glob('*.npz')]}
 for row in a['states']:
  s=lookup[row['step']];g=corotated(r,s['x'],s['d'])
  for node in row['local_nodes']:
   q=node['atom_id']-1;condition=float(g['condition'][q]);node['local_affine_fit_condition']=condition;node['local_affine_detF']=float(g['detF'][q]) if condition<1e10 else None;node['local_affine_fit_status']='DIAGNOSTIC ONLY' if condition<1e10 else 'NOT MEASURED'
 note=' Local det(F) comes from a damaged-neighborhood affine least-squares fit; with strong non-affinity it is not a measured continuum Jacobian or material density.'
 if note not in a['local_density_limit']:a['local_density_limit']+=note
 dump(p,a)
 p=O/'MIXED_MODE_DIAGNOSTIC.json';a=json.loads(p.read_text())
 for frame in a['frames']:
  for b in frame['candidate_bridges']:
   b['shear_status']='DIAGNOSTIC ONLY' if b['valid'] else 'NOT MEASURED'
   if not b['valid'] and b['shear'] is not None:b['invalid_raw_shear_not_for_interpretation']=b['shear'];b['shear']=None
 dump(p,a)
 p=O/'V3_CHUNK_CANDIDATES.json';a=json.loads(p.read_text());a['dominant_attached_body_screen']=[]
 for path,s in keyframes():
  lab,sz,att=topology(r,s['d']);q=np.flatnonzero(att)[np.argmax(sz[att])];nodes=np.flatnonzero(lab==q);out,inside=boundary(r,nodes);inside=inside[r['bond_types'][inside]==1]
  a['dominant_attached_body_screen'].append(dict(step=int(s['step']),nodes=len(nodes),external_reference_bonds=len(out),external_alive_bulk=int(np.sum((r['bond_types'][out]==1)&(s['d'][out]<1))),external_alive_interface=int(np.sum((r['bond_types'][out]==2)&(s['d'][out]<1))),mean_external_d=float(np.mean(s['d'][out])),mean_internal_d=float(np.mean(s['d'][inside])),internal_survival=float(np.mean(s['d'][inside]<1))))
 a['search_scope_detail']='Local candidate proposals cover 3..963 material points. The entire dominant attached body is additionally tabulated; arbitrary >963-node subregions were not exhaustively enumerated.';dump(p,a)
 gpath=R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/V3_MATERIAL_POINTS_AND_BONDS.npz';g=load_npz(gpath)
 for k in ['X','bonds','bond_types']:assert np.array_equal(g[k],r[k])
 dump(O/'EXPOSED_SURFACE_SOURCE_IDENTITY.json',dict(status='PASS',geometry=str(gpath.relative_to(R)),geometry_sha256=sha(gpath),reference_X_bonds_and_types_exact=True,builder='pd_clot/lammps_fragmentation/y_clot_morphology_v3.py::geometry -> y_mural_clot_geometry.py::build -> exposed_surface.vtp'))
 print('AUGMENTED_AUDITS')
if __name__=='__main__':main()
