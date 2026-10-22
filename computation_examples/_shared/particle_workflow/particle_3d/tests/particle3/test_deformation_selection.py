import json
import numpy as np
import pytest
from particle_3d.particle3_cases import cylinder_wall
from particle_3d.rbc_capillary_surrogate import select_rbc_shape,area_feasible_interval
from particle_3d.particle_shapes import Capsule,Ellipsoid,EPS


@pytest.mark.parametrize('index',range(5))
def test_trigger_minimal_deformation_and_explicit_infeasible(p3_geometries,index):
    g=p3_geometries[index];rmin,req,budget=area_feasible_interval(g)
    for radius in [g.a_m*1.2,(rmin+min(req,g.a_m))*.5,.8*rmin]:
        wall=cylinder_wall(radius);decision=select_rbc_shape([0,0,0],g,[1,0,0,0],[0,0,1e-4],wall)
        apothem=wall.provenance['tube_apothem_m']
        if apothem>=g.a_m:
            assert decision.status=='FREE_OBLATE' and isinstance(decision.shape,Ellipsoid)
        elif apothem<rmin:
            assert decision.status=='DEFORMATION_SURROGATE_INFEASIBLE' and decision.shape is None
        else:
            assert decision.status=='CAPILLARY_DEFORMED' and decision.original_oblate_gap_m<0
            c=decision.shape
            assert abs(c.radius_m-min(req,apothem))<1024*EPS*req
            assert abs(c.volume_m3/g.volume_m3-1)<16*EPS
            assert c.area_m2<=budget*(1+512*EPS)
            assert decision.gap.gap_m>=-decision.gap.roundoff_m
            np.testing.assert_array_equal(c.axis_world,[0,0,1])


def test_zero_flow_and_search_guard_are_explicit(p3_geometries,real_wall,p3_repo):
    g=p3_geometries[0];wall=cylinder_wall(1.5e-6)
    assert select_rbc_shape([0,0,0],g,[1,0,0,0],[0,0,0],wall).status=='DEFORMATION_AXIS_UNRESOLVED'
    initial=json.loads((p3_repo/'particle_3d/reports/particle1/data/05_real_initialization.json').read_text())
    from particle_3d.particle2_cases import initial_orientation
    decision=select_rbc_shape(initial['initial_position_m'],g,initial_orientation(),initial['initial_velocity_m_s'],real_wall,node_limit=1)
    assert decision.status=='DEFORMATION_SEARCH_UNRESOLVED' and decision.shape is None


def test_all_1152_distribution_cases_conserve_volume_area_and_radius(p3_repo):
    rows=json.loads((p3_repo/'particle_3d/reports/particle3/data/07_feasibility_map.json').read_text())
    assert len(rows)==128*9
    assert {r['status'] for r in rows}=={'FREE_OBLATE','CAPILLARY_DEFORMED','DEFORMATION_SURROGATE_INFEASIBLE'}
    for row in rows:
        if row['status']=='CAPILLARY_DEFORMED':
            assert row['volume_relative_error']<16*EPS
            assert row['area_violation_m2']<=512*EPS*row['area_budget_m2']
            assert abs(row['R_cap_m']-row['analytic_max_radius_m'])<1024*EPS*row['analytic_max_radius_m']


def test_quasistatic_path_proof_keeps_actual_maximal_radius(p3_repo,p3_geometries,real_wall):
    from particle_3d.field import FrozenFEMField
    from particle_3d.particle2_cases import initial_orientation
    from particle_3d.particle3_motion import initial_wall_state,contact_trial
    from particle_3d.physical_time_refinement import refine_interval,swept_clearance_certificate
    field=FrozenFEMField.from_frozen(p3_repo/'formal_3D_flow_solver/FEM_SimVascular')
    initial=json.loads((p3_repo/'particle_3d/reports/particle1/data/05_real_initialization.json').read_text())
    g=p3_geometries[0];start=initial_wall_state(initial['initial_position_m'],initial_orientation(),field,real_wall,geometry=g)
    end,ledger=refine_interval(start,initial['validation_timesteps_s'][-1],contact_trial(field,real_wall))
    proof=ledger[0]['proof'];assert 'PROOF_ONLY_WITNESS_R_M=' in proof
    radius=float(proof.split('PROOF_ONLY_WITNESS_R_M=')[1]);rmin,_,budget=area_feasible_interval(g)
    assert rmin<=radius<min(start.shape.radius_m,end.shape.radius_m)
    witness=[]
    for state in [start,end]:
        c=Capsule(state.shape.center_m,state.shape.axis_world,radius,g.volume_m3/(np.pi*radius**2)-4*radius/3)
        assert c.area_m2<=budget*(1+512*EPS)
        assert abs(c.volume_m3/g.volume_m3-1)<16*EPS
        witness.append(c)
    assert swept_clearance_certificate(*witness,real_wall)[0]
    selected=select_rbc_shape(end.shape.center_m,g,end.quaternion_wxyz,end.free_velocity_m_s,real_wall)
    assert selected.shape.radius_m==end.shape.radius_m
    np.testing.assert_array_equal(start.quaternion_wxyz,end.quaternion_wxyz)
    np.testing.assert_allclose(end.shape.axis_world,end.free_velocity_m_s/np.linalg.norm(end.free_velocity_m_s),rtol=0,atol=2*EPS)


@pytest.mark.parametrize('index,center',[
    (0,[0.00010339739355737996,5.519822873878984e-05,0.00015579560614967613]),
    (1,[0.00010332635655586305,5.3265641582302286e-05,0.00015366561737568528])])
def test_search_refines_small_interval_above_best_until_global_bound_closes(p3_repo,p3_geometries,real_wall,index,center):
    from particle_3d.field import FrozenFEMField
    from particle_3d.particle2_cases import initial_orientation
    from particle_3d.particle_shapes import roundoff_length
    field=FrozenFEMField.from_frozen(p3_repo/'formal_3D_flow_solver/FEM_SimVascular')
    g=p3_geometries[index]
    decision=select_rbc_shape(center,g,initial_orientation(),field.sample(center).velocity_m_s,real_wall)
    assert decision.status=='CAPILLARY_DEFORMED',decision.reason
    assert decision.radius_optimality_bound_m<=roundoff_length(area_feasible_interval(g)[1])
    assert decision.gap.gap_m>=-decision.gap.roundoff_m


def test_smaller_radius_can_fail_by_lengthening_in_finite_geometry(p3_geometries):
    from particle_3d.particle3_cases import plane_triangle
    from particle_3d.rbc_orientation import rotation_matrix,quaternion_from_short_axis
    from particle_3d.wall_geometry import WallGeometry
    from particle_3d.wall_gap import wall_gap
    from particle_3d.rbc_capillary_surrogate import capsule_length
    g=p3_geometries[0];side=1.6e-6
    end=1.5e-6+capsule_length(g.volume_m3,1.5e-6)/2
    triangles=[]
    for direction,offset in [([1,0,0],-side),([-1,0,0],-side),([0,0,1],-end),([0,0,-1],-end)]:
        n=np.array(direction);rotation=rotation_matrix(quaternion_from_short_axis(n))
        triangles.append(plane_triangle()@rotation.T+offset*n)
    wall=WallGeometry(triangles)
    def capsule(radius):return Capsule([0,0,0],[0,0,1],radius,capsule_length(g.volume_m3,radius))
    assert wall_gap(capsule(1.1e-6),wall).state=='PENETRATING'  # thinner but too long
    assert wall_gap(capsule(1.55e-6),wall).state=='SEPARATED'
    assert wall_gap(capsule(1.7e-6),wall).state=='PENETRATING'  # too wide
    d=select_rbc_shape([0,0,0],g,[1,0,0,0],[0,0,1e-4],wall)
    assert d.status=='CAPILLARY_DEFORMED'
    assert abs(d.shape.radius_m-side)<=d.radius_optimality_bound_m


@pytest.mark.parametrize('index',range(2))
def test_cached_radius_certificates_match_uncached_exact_search(p3_repo,p3_geometries,real_wall,index):
    from particle_3d.field import FrozenFEMField
    from particle_3d.particle2_cases import initial_orientation
    field=FrozenFEMField.from_frozen(p3_repo/'formal_3D_flow_solver/FEM_SimVascular')
    positions=[[0.00010343884423491545,5.7755787565838546e-05,0.00015849126066314057],
               [0.00010339739355737996,5.519822873878984e-05,0.00015579560614967613]]
    g=p3_geometries[index]
    for center in positions:
        args=(center,g,initial_orientation(),field.sample(center).velocity_m_s,real_wall)
        cached=select_rbc_shape(*args);reference=select_rbc_shape(*args,use_radius_cache=False)
        assert cached.status==reference.status=='CAPILLARY_DEFORMED'
        assert abs(cached.shape.radius_m-reference.shape.radius_m)<=cached.radius_optimality_bound_m+reference.radius_optimality_bound_m
        assert abs(cached.gap.gap_m-reference.gap.gap_m)<=cached.gap.roundoff_m
