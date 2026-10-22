"""Small actual-field subsets for independent WSS review without whole CFD files."""
import argparse
import numpy as np
import pyvista as pv
from case_common import *
from flow_solver_support.wss_case import material
from vessel_regions import CENTERS
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);a=p.parse_args();case=a.case.resolve()
assert json.loads((case/'reports/flow_quality.json').read_text())['accepted_final_and_log_checks']
m=np.load(case/'SV_MESH/mesh_arrays.npz');f=np.load(case/'frozen_flow/flow_arrays_si.npz')
x=m['points_m'];t=m['tetra'];b=m['boundary_triangles'];tags=m['facet_tags']
surface=pv.read(case/'wss/data/wall_wss_si.vtp');ids=np.asarray(surface.cell_data['Global_boundary_facet_zero_based'])
cent=x[b[ids]].mean(axis=1);selected=np.logical_or.reduce([np.linalg.norm(cent-np.array(CENTERS[n])*1e-6,axis=1)<5e-6 for n in ['J1_r5um','J2_r5um']])
facets=ids[selected];owners=np.asarray(surface.cell_data['Parent_tetra_zero_based'])[selected];cells=np.unique(owners);nodes,local_t=np.unique(t[cells],return_inverse=True)
local_b=np.searchsorted(nodes,b[facets]);assert np.array_equal(nodes[local_b],b[facets])
out=V/'evidence/local_wall_snapshots';out.mkdir(parents=True,exist_ok=True)
np.savez_compressed(out/(case.name+'.npz'),points_m=x[nodes],tetra=local_t.reshape(-1,4),wall_triangles=local_b,velocity_m_s=f['velocity_m_s'][nodes],pressure_Pa=f['pressure_pa'][nodes],original_node_zero_based=nodes,original_tetra_zero_based=cells,original_wall_facet_zero_based=facets,expected_WSS_raw_Pa=np.asarray(surface.cell_data['WSS_raw_Pa'])[selected],expected_traction_Pa=np.asarray(surface.cell_data['Tangential_viscous_traction_Pa'])[selected],mu_Pa_s=np.array(material(case)['mu_Pa_s']))
dump(out/(case.name+'.json'),dict(case=case.name,source_case=str(case),core_module_sha256=sha(wss.__file__),source_mesh_sha256=sha(case/'SV_MESH/mesh_arrays.npz'),source_flow_sha256=sha(case/'frozen_flow/flow_arrays_si.npz'),source_wall_sha256=sha(case/'wss/data/wall_wss_si.vtp'),subset_sha256=sha(out/(case.name+'.npz')),selection='Union of J1 and J2 radius5um wall-facet-centroid masks',nodes=len(nodes),tetra=len(cells),wall_facets=len(facets),description='Actual unchanged nodal CFD values and original coordinates; subset is NOT a closed CFD domain and is NOT a new solve. Used only to independently repeat local production WSS recovery.'))
print(case.name,len(nodes),'nodes',len(facets),'wall facets', (out/(case.name+'.npz')).stat().st_size,'bytes')
