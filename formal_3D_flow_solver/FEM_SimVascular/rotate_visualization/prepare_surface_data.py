"""Replace derived surface data using the existing P1 WSS implementation."""
from pathlib import Path
import importlib.util,json,sys,time
import numpy as np
import pyvista as pv

HERE=Path(__file__).resolve().parent
CASE=HERE/'input_data'

def main():
    started=time.perf_counter()
    context=json.loads((HERE/'BUILD_CONTEXT.json').read_text())
    source=Path(context['source_case']);flow=source.parents[1]
    implementation=flow/'scripts/flow_2mmps/compute_field_diagnostics.py'
    spec=importlib.util.spec_from_file_location('existing_field_diagnostics',implementation)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    sys.path.insert(0,str(flow/'src'))
    from sv_validation.postprocess import SolutionMeasurements
    manifest=json.loads((source/'frozen_flow/manifest.json').read_text())
    for name,entry in manifest['files'].items():
        assert old.sha(source/'frozen_flow'/name)==entry['sha256'],name
    assert old.sha(CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')==context['source_field_sha256']
    assert old.sha(CASE/'frozen_flow/flow_arrays_si.npz')==old.sha(source/'frozen_flow/flow_arrays_si.npz')
    policy=json.loads((source/'policy.json').read_text())
    measurements=SolutionMeasurements(source/'SV_MESH/mesh_arrays.npz',policy['Q_target_m3_s'],policy['Umean_m_s'])
    u,p=measurements.read(CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    actual=measurements.measure(u,p)
    assert actual['epsilon_Q']<=1e-6 and actual['epsilon_mass']<=1e-6
    assert actual['wall_noslip_pass'] and actual['velocity_finite'] and actual['pressure_finite']
    frozen=np.load(CASE/'frozen_flow/flow_arrays_si.npz')
    points,tetra,triangles,tags,velocity,pressure=[frozen[k] for k in ['points_m','tetra','boundary_triangles','facet_tags','velocity_m_s','pressure_pa']]
    assert np.array_equal(u,velocity) and np.array_equal(p,pressure)
    assert np.array_equal(points,measurements.points) and np.array_equal(tetra,measurements.tetra)
    wall_ids=np.flatnonzero(tags==1);wall=triangles[wall_ids]
    owners=old.boundary_owners(tetra,wall)
    centers,area,normal=old.wall_geometry(points,tetra,wall,owners)
    gradient=old.p1_gradients(points,tetra,velocity)
    viscosity=.00345312
    traction=old.tangential_traction(gradient[owners],normal,viscosity)
    magnitude=np.linalg.norm(traction,axis=1)
    display,weight=old.nodal_average(wall,magnitude,area,len(points))
    mesh,ids=old.surface(points,wall)
    mesh.point_data['WSS_display_Pa']=display[ids]
    mesh.point_data['Pressure_Pa']=pressure[ids]
    for key,value in {'WSS_raw_Pa':magnitude,'Tangential_viscous_traction_Pa':traction,
                      'Outward_normal':normal,'Area_m2':area,'Parent_tetra_zero_based':owners,
                      'Global_boundary_facet_zero_based':wall_ids}.items():mesh.cell_data[key]=value
    out=CASE/'field_diagnostics';mesh.save(out/'data/wall_wss_si.vtp')
    full,pressure_ids=old.surface(points,triangles)
    full.point_data['Pressure_Pa']=pressure[pressure_ids];full.cell_data['Boundary_tag']=tags
    full.save(out/'data/pressure_surface_si.vtp')
    selected=np.unique(owners)[np.linspace(0,len(np.unique(owners))-1,128,dtype=int)]
    affine=np.concatenate([np.ones((len(selected),4,1)),points[tetra[selected]]],axis=2)
    independent=np.swapaxes(np.linalg.solve(affine,velocity[tetra[selected]])[:,1:],1,2)
    gradient_error=float(abs(independent-gradient[selected]).max())
    assert np.allclose(independent,gradient[selected],rtol=2e-10,atol=1e-7)
    tangent_error=float(abs(np.einsum('ij,ij->i',traction,normal)).max())
    assert tangent_error<max(1.,magnitude.max())*1e-12
    assert np.isfinite(magnitude).all() and np.all(weight[ids]>0)
    assert np.linalg.norm(velocity[ids],axis=1).max()==0
    record=dict(all_pass=True,case=source.name,source_field_sha256=context['source_field_sha256'],
        source_case=str(source),original_diagnostics_implementation=str(implementation),
        original_diagnostics_sha256=old.sha(implementation),viscosity_Pa_s=viscosity,
        pressure_Pa=dict(min=float(pressure.min()),max=float(pressure.max()),surface_min=float(pressure[pressure_ids].min()),surface_max=float(pressure[pressure_ids].max())),
        wss_Pa=dict(raw_min=float(magnitude.min()),raw_max=float(magnitude.max()),
            area_weighted_mean=float(np.sum(area*magnitude)/area.sum()),display_min=float(display[ids].min()),display_max=float(display[ids].max())),
        wss_method='Exact P1 tetra velocity gradients; tangential viscous traction on outward planar WALL facets',
        wss_display='Area-weighted mean of incident facet magnitudes at original wall nodes; geometry unchanged',
        wss_vector_convention='(I-n*n^T)*mu*(grad(u)+grad(u)^T)*n',
        wss_not_native_saved_solver_output=True,mesh_convergence_not_assessed=True,
        wall_facets=len(wall),independent_gradient_max_absolute_error_s_inv=gradient_error,
        tangency_max_error_Pa=tangent_error,unique_parent_for_all_wall_facets=True,max_wall_speed_m_s=0.,
        physical_geometry_unchanged=True,source_values_clipped=False,
        outputs_sha256={str(f.relative_to(out)):old.sha(f) for f in (out/'data').glob('*.vtp')},
        elapsed_seconds=time.perf_counter()-started)
    old.dump(out/'COMPUTE_VALIDATION.json',record)
    old.dump(HERE/'audit/NEW_FLOW_REVIEW.json',dict(all_pass=True,source_field_sha256=context['source_field_sha256'],
        frozen_manifest_hashes_match=True,NPZ_and_VTU_values_exactly_equal=True,measurements=actual,
        pressure=record['pressure_Pa'],WSS=record['wss_Pa'],source_solver_case_unchanged=True))
    print(json.dumps(record,indent=2),flush=True)

if __name__=='__main__':main()
