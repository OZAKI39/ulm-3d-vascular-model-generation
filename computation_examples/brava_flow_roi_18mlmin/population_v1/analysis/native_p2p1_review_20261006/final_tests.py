"""Final engineering/regression checks; scientific geometry gate remains blocked."""
from common import *
from p2_sampler import *
from PIL import Image

def main():
    tests=list(read(HERE/'P2_FIELD_SAMPLER_TESTS.json')['tests'])
    def check(name,condition,**details):tests.append(dict(name=name,status='PASS' if condition else 'FAIL',**details))
    b=read(HERE/'P2_REPRESENTATION_FIXTURES.json');d=read(HERE/'P2_POISEUILLE_REPRESENTATION.json');g=read(HERE/'P2_CURVED_GEOMETRY_AUDIT.json')
    unit=read(HERE/'NATIVE_UNIT_FIXTURE.json')
    check('native_mixed_assembly_dimensions',unit['local_residual_shape']==[4,10] and unit['local_tangent_shape']==[16,10,10] and unit['mixed_edge_coupling_abs_sum']>0)
    check('native_pressure_edge_identity_rows',unit['native_pressure_edge_handling_pass'] and unit['edge_pressure_residual_abs_sum']==0)
    check('native_linear_continuity_residual',unit['continuity_linear_field_error']<1e-11)
    check('complete_B1_quadratic_reproduction',b['B1']['status']=='PASS')
    check('B2_wall_edge_BC_and_rank',b['B2']['status']=='PASS' and len(b['B2']['wall_local_nodes'])==6 and b['B2']['P1_gradient_ranks']==[1] and b['B2']['P2_gradient_ranks']==[3])
    target=b['real_138584'];edges=target['edges']
    check('real_138584_interior_edge_ownership',len(edges)==6 and sum(r['fixed_velocity'] for r in edges)==5 and target['free_P2_velocity_nodes']==1)
    check('real_138584_no_invented_internal_DOF',target['unit_free_DOF_fixture']['nonzero_internal_samples']==4096 and target['unit_free_DOF_fixture']['gradient_ranks']==[1])
    check('D1_all_frozen_and_random_samples',d['frozen_probes']==960 and d['random_points']==4096 and d['invalid']==0)
    check('D1_roundoff_Poiseuille',d['D1']=='PASS' and max(d['physical_scale_max_errors'].values())<1e-11)
    field=np.load(HERE/'runtime/PIPE_AFFINE_ELEVATION.npz')
    raw=np.load(BASE/'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu/SV_MESH/mesh_arrays.npz')
    X,C,edges=field['points_m'],field['tetra10'],field['unique_edges'];n=len(raw['points_m'])
    check('elevation_corner_positions_unchanged',np.array_equal(X[:n],raw['points_m']))
    check('elevation_original_cell_topology_unchanged',np.array_equal(np.sort(C[:,:4],axis=1),np.sort(raw['tetra'],axis=1)))
    check('elevation_unique_edge_numbering',len(np.unique(edges,axis=0))==len(edges) and np.array_equal(X[n:],X[edges].mean(1)))
    pairs=np.sort(C[:,:4][:,EDGES],axis=2)
    check('global_shared_midpoint_consistency',np.array_equal(edges[C[:,4:]-n],pairs))
    check('boundary_tag_inheritance',np.array_equal(field['facet_tags'],raw['facet_tags']) and np.array_equal(field['boundary6'][:,:3],raw['boundary_triangles']))
    B=field['boundary6'];epairs=np.sort(B[:,:3][:,[[0,1],[1,2],[2,0]]],axis=2)
    check('boundary_midpoint_connectivity',np.array_equal(edges[B[:,3:]-n],epairs))
    check('affine_positive_Jacobians',d['mesh']['min_positive_Jacobian']>0)
    for f in g['fixtures']:
        check(f['label']+'_positive_curved_Jacobian_all_reference_points',f['all_reference_Jacobians_positive_certificate']['passed'])
        check(f['label']+'_wall_and_cap_projection',f['wall_midpoint_radial_error_m']<1e-14 and (f['inlet_midpoint_plane_error_m'] or 0)<1e-14)
        check(f['label']+'_affine_control_roundoff',f['variants']['affine']['continuity_geometry_weight_relative_L2_mismatch']<1e-11)
        check(f['label']+'_native_construct_fluid_matches_audited_control_flow',f['variants']['curved']['native_prediction_relative_max_error']<1e-11)
        ratios=f['variants']['curved']['continuity_used_over_velocity_Jacobian']
        check(f['label']+'_detect_curved_weight_inconsistency',abs(ratios[0]-1)>1e-6 and np.max(abs(np.array(ratios[1:])-1))<1e-11,
              meaning='Test PASS confirms a scientific implementation blocker; does not certify CFD.')
    baseline=read(HERE/'INPUT_SHA256.json')
    mismatch=[k for k,h in baseline.items() if sha(BASE/k)!=h]
    check('historical_inputs_P1_fields_and_contracts_unchanged',not mismatch,mismatches=mismatch,files=len(baseline))
    previous=read(PREVIOUS/'EXACT_DEPLOYED_SOURCE_MANIFEST.json')
    check('production_retained_source_unchanged',all(sha(RETAINED/k)==h for k,h in previous.items()),files=len(previous))
    source=read(HERE/'SOURCE_IDENTITY.json')
    check('production_binary_and_PETSc_unchanged',sha(source['binary'])==source['binary_sha256'] and sha(source['PETSc']['path'])==source['PETSc']['sha256'])
    for p in sorted(HERE.glob('0*.png')):
        with Image.open(p) as im:check(p.stem+'_300dpi',abs(im.info['dpi'][0]-300)<.1)
    check('gate_blocks_CFD_without_alternate_route',g['status']=='P2_CURVED_GEOMETRY_NOT_SUPPORTED' and not (HERE/'pipe_CFD').exists() and not (HERE/'BraVa_CFD').exists())
    result=dict(status='PASS' if all(t['status']=='PASS' for t in tests) else 'FAIL',passed=sum(t['status']=='PASS' for t in tests),
        failed=sum(t['status']=='FAIL' for t in tests),tests=tests,
        scope='Engineering correctness and detection of curved implementation inconsistency; not CFD physical accuracy.',
        not_run=['GPU P2 CFD','P2 production particle sampler','trajectory replay','full curved mesh generation'])
    put('TEST_RESULTS.json',result)
    print(result['status'],result['passed'],'PASS',result['failed'],'FAIL')
    if result['failed']:raise SystemExit(1)

if __name__=='__main__':main()
