"""Read-only analysis of existing tracks; static terminal solves never advance time."""
from pathlib import Path
import collections,csv,gzip,hashlib,json,os,sys,time
import numpy as np
from scipy.optimize import linprog,nnls
import pyvista as pv

OUT=Path(__file__).resolve().parents[1];RUN=OUT.parent;PROJECT=RUN.parents[2]
sys.path.insert(0,str(PROJECT/'particle_3d/src'))
from particle_3d.field import FrozenFEMField
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_gap import wall_gap
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle9a_motion import augment_planar_system
from particle_3d.resistance_solver import solve_resistance,contact_jacobian
from particle_3d.open_boundary_rim import build_rim_topology

VIS=Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1')
PATHS=VIS.parents[1]/'wss_validation_v2/meeting_question_maps/data/reference_paths_um.csv'
J1=np.array([92.,49.,111.]);J2=np.array([130.04,82.04,87.18])
LOCK={}

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def read(p):
    p=Path(p);LOCK[str(p)]=sha(p);return json.loads(p.read_text())

def dump(name,v):
    (OUT/name).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False,
        default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())+'\n')

def csvout(name,rows):
    with (OUT/'data'/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def distribution(v):
    v=np.asarray(v,float)
    return dict(min=float(v.min()),median=float(np.median(v)),max=float(v.max()),mean=float(v.mean()))

def branch_at(point,paths):
    best=None
    for name,rr in paths.items():
        xyz=np.array([[float(r[k]) for k in ['x_um','y_um','z_um']] for r in rr])
        arc=np.array([float(r['s_um']) for r in rr]);v=np.diff(xyz,axis=0);length=np.linalg.norm(v,axis=1)
        good=length>1e-12;xyz0=xyz[:-1][good];v=v[good];length=length[good]
        t=np.clip(np.einsum('ij,ij->i',point-xyz0,v)/length**2,0,1)
        d=np.linalg.norm(point-xyz0-t[:,None]*v,axis=1);k=int(d.argmin())
        item=(float(d[k]),name,float(arc[:-1][good][k]+t[k]*length[k]))
        if best is None or item<best:best=item
    return best

def lp_direction(normals,u):
    unit=u/np.linalg.norm(u)
    lp=linprog(-unit,A_ub=-normals,b_ub=np.zeros(len(normals)),bounds=[(-1,1)]*3,method='highs',
        options=dict(primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10))
    assert lp.success
    lam,res=nnls(normals.T,-unit)
    value=float(unit@lp.x);tol=1e-8
    if value>tol:state='DOWNSTREAM_DIRECTION_EXISTS'
    elif res<tol:state='NO_DOWNSTREAM_DIRECTION_CERTIFIED'
    else:state='UNRESOLVED'
    # Nonzero escape/retraction directions may exist even if forward is blocked.
    axes=[]
    for axis in np.eye(3):
        for sign in [-1,1]:
            t=linprog(-sign*axis,A_ub=-normals,b_ub=np.zeros(len(normals)),bounds=[(-1,1)]*3,method='highs')
            assert t.success;axes.append(float(sign*axis@t.x))
    return dict(status=state,maximum_forward_projection=value,dual_certificate_residual=float(res),
        direction=lp.x,minimum_normal_projection=float((normals@lp.x).min()),
        any_nonzero_feasible_direction=max(axes)>tol,tolerance=tol)

def main():
    start=time.time()
    summary=read(RUN/'data/final_summary.json');rows=read(RUN/'data/metrics.json');cohort=read(RUN/'data/cohort.json')
    events={e['particle_id']:e for e in cohort['events']};assert len(rows)==1500
    assert summary['flow_sha256']=='fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12'
    sources=summary['identity']['protected_sources']
    for relative,h in sources.items():assert sha(PROJECT/relative)==h,relative
    for relative,h in summary['identity']['acceleration_sources'].items():assert sha(RUN/'scripts'/relative)==h,relative
    for f in ['formal_dynamics_p9a5.py','formal_cohort_p9a5.py','particle82a_integration.py','particle81_simulation.py',
              'nearfield_regularization.py','nearfield_handoff.py','resistance_solver.py','particle9a_motion.py',
              'planar_wall_hydrodynamics.py','field.py','convex_triangle.py','wall_gap.py']:
        p=PROJECT/'particle_3d/src/particle_3d'/f;LOCK[str(p)]=sha(p)
    for name in ['scripts/render_results.py','additional_views/render_views.py','scripts/campaign.py','data/render_manifest.json']:
        LOCK[str(RUN/name)]=sha(RUN/name)
    delivery=read(RUN/'DELIVERY_COMPLETE.json');meshfile=RUN/'data/gpu_mesh_input.npz'
    assert sha(meshfile)==delivery['files']['data/gpu_mesh_input.npz'];LOCK[str(meshfile)]=sha(meshfile)
    arrays=np.load(meshfile)
    flowfile=VIS/'input_data/frozen_flow/steady_flow_mean_2p0_mmps.vtu';assert sha(flowfile)==summary['flow_sha256'];LOCK[str(flowfile)]=sha(flowfile)
    flow=pv.read(flowfile)
    assert np.array_equal(flow.points,arrays['points']) and np.array_equal(flow['Velocity'],arrays['velocity'])
    assert np.array_equal(flow['Pressure'],arrays['pressure'])
    field=FrozenFEMField(arrays['points'],arrays['tetra'],arrays['velocity'],arrays['pressure'])
    wall=WallGeometry(arrays['wall_triangles'],provenance=dict(role='EXACT_CURRENT_PRODUCTION_WALL_ARRAY_SNAPSHOT',sha256=sha(meshfile)))
    wallfile=VIS/'input_data/solver_mesh/mesh-surfaces/WALL.vtp';wm=pv.read(wallfile);LOCK[str(wallfile)]=sha(wallfile)
    faces=wm.faces.reshape(-1,4)[:,1:]
    assert np.array_equal(wm.points[faces],wall.triangles)
    wall.global_node_ids=wm.point_data['GlobalNodeID'][faces]-1
    caps={};cap_centers={}
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        p=wallfile.parent/f'{role}.vtp';cap=pv.read(p);LOCK[str(p)]=sha(p)
        cf=cap.faces.reshape(-1,4)[:,1:];caps[role]=cap.point_data['GlobalNodeID'][cf]-1
        cap_centers[role]=np.array(cap.center)*1e6
    wall.open_boundary_topology=build_rim_topology(wall.global_node_ids,caps)
    LOCK[str(PATHS)]=sha(PATHS);pathrows=list(csv.DictReader(PATHS.open()));paths={}
    for r in pathrows:paths.setdefault(r['path'],[]).append(r)
    policy=NearFieldRegularizationV1();mu=.00345312
    allrows=[];stops=[];contacts_csv=[];static=[];representative={};verified_files=0
    safety=collections.Counter();same_last_windows=0;stopped_raw={}
    for index,r in enumerate(rows):
        pid=r['particle_id'];folder=RUN/'tracks'/f'mb_{pid:06d}'
        marker=read(folder/'COMPLETE.json');assert marker['identity']==summary['identity']
        # Hash all original track files, including compressed audit streams.
        for name,h in marker['files'].items():
            path=folder/name;actual=sha(path);assert actual==h,str(path);LOCK[str(path)]=actual;verified_files+=1
        metric=read(folder/'metrics.json');meta=read(folder/'trajectory.json')
        for key in ['status','outlet','end_reason','failure_detail','diameter_um','final_position_m']:assert metric[key]==r[key]
        samples=np.load(folder/'trajectory.npz')['samples'];assert np.isfinite(samples).all()
        assert np.all(np.diff(samples[:,0])>0) and np.diff(samples[:,0]).max()<=.0005+1e-12
        assert meta['formal_identity']==summary['identity'] and meta['integration_config']['dt_s']==.0005
        assert np.array_equal(samples[0,1:4],events[pid]['birth_center_m'])
        tail=samples[-1];stationary=r['status']=='SUPPORTED_STATIONARY'
        after=collections.Counter({k:r[k] for k in ['penetration_count','handoff_violation_count','inlet_escape_count','nan_inf_count','unclassified_corruption_count']});safety.update(after)
        base=dict(particle_id=pid,status=r['status'],outlet=r['outlet'] or '',diameter_um=r['diameter_um'],
            termination_age_ms=tail[0]*1e3,final_x_um=tail[1]*1e6,final_y_um=tail[2]*1e6,final_z_um=tail[3]*1e6,
            final_particle_speed_mm_s=float(np.linalg.norm(tail[4:7])*1e3),
            path_length_um=r['path_length_m']*1e6,wall_gap_final_nm=float(tail[14]*1e9),
            failure_detail=r['failure_detail'] or '',sample_count=len(samples))
        allrows.append(base)
        if not stationary:continue
        stopped_raw[pid]=samples
        support=read(folder/'support.json');window=support['last_32_accepted_contact_solves'];assert len(window)==32
        assert all(a['solver']['contact_redundancy']['rank_after']==3 and np.max(abs(np.array(a['velocity'][:3])))<=a['solver']['contact_kkt']['velocity_budget_m_s']
                   and np.min(a['solver']['multipliers'])>=0 for a in window)
        # Independently stream the saved audit and compare its final 32 solves.
        last=collections.deque(maxlen=32)
        with gzip.open(folder/'audit.jsonl.gz','rt') as f:
            for line in f:
                a=json.loads(line)
                if a['accepted'] and 'solver' in a:last.append(a)
        assert list(last)==window;same_last_windows+=1
        final=window[-1];position=np.array(final['position']);shape=Sphere(position,r['radius_m']);shapes={pid:shape}
        sample=field.sample(position);assert sample.inside_lumen
        free=np.r_[sample.velocity_m_s,.5*sample.vorticity_s_inv]
        assert np.max(abs(free-np.array(final['free'])))<1e-8
        current_contacts,records=handoff_constraints(shapes,wall,policy)
        system=assemble_v1(shapes,{pid:free},mu,wall,policy=policy)
        system=augment_planar_system(system,shapes,mu,wall,lambda p:field.sample(p).velocity_gradient_s_inv,wall_gap)
        solution=solve_resistance(system,particles=shapes,constraints=current_contacts)
        j,ordered=contact_jacobian(system,shapes,current_contacts);normals=j[:,:3]
        assert [list(a.canonical_id) for a in ordered]==final['solver']['contact_ids']
        verr=float(np.max(abs(solution.velocity[:3]-np.array(final['velocity'][:3]))))
        herr=float(np.max(abs(solution.unconstrained[:3]-np.array(final['hydro'][:3]))))
        assert verr<1e-12 and herr<1e-12
        lp=lp_direction(normals,sample.velocity_m_s)
        plane=wall_gap(shape,wall);lower=policy.lower_handoff_gap(shape.radius_m)['h_lower_m']
        assert abs(plane.gap_m-lower)<=wall.roundoff_m
        for cc,record in zip(current_contacts,records):
            if record['contact_type']!='CONTINUUM_HANDOFF_CONTACT':raise ValueError('Unexpected exact hard contact')
            contacts_csv.append(dict(particle_id=pid,triangle_id=cc.canonical_id[3],feature=cc.canonical_id[4],
                h_geom_nm=record['h_geom_m']*1e9,h_lower_nm=record['h_lower_m']*1e9,
                nx=cc.normal[0],ny=cc.normal[1],nz=cc.normal[2],
                wall_x_um=record['wall_point_m'][0]*1e6,wall_y_um=record['wall_point_m'][1]*1e6,wall_z_um=record['wall_point_m'][2]*1e6))
        last64=samples[-64:];span=float(np.linalg.norm(np.ptp(last64[:,1:4],axis=0)))
        nominal=np.arange(0,tail[0]+1e-14,.0005)
        idx=np.array([np.argmin(abs(samples[:,0]-t)) for t in nominal]);nt=samples[idx]
        equal=np.all(np.diff(nt[:,1:4],axis=0)==0,axis=1);trailing=0
        for b in equal[::-1]:
            if not b:break
            trailing+=1
        if '32_NOMINAL' in r['failure_detail']:assert trailing>=32
        else:assert span<=16*wall.roundoff_m
        near,path,s=branch_at(tail[1:4]*1e6,paths)
        distj1=float(np.linalg.norm(tail[1:4]*1e6-J1));distj2=float(np.linalg.norm(tail[1:4]*1e6-J2))
        region='J1_r5um' if distj1<5 else 'J2_r5um' if distj2<5 else path
        # Onset below 1e-9 m/s is a retrospective diagnostic, not a solver stop rule.
        speed=np.linalg.norm(samples[:,4:7],axis=1);moving=np.flatnonzero(speed>1e-9)
        onset=float(samples[moving[-1],0]) if len(moving) else 0.
        cosine=normals@normals.T;angle=float(np.degrees(np.arccos(np.clip(cosine,-1,1))).max())
        extra=dict(region=region,nearest_path=path,path_arc_um=s,distance_to_path_um=near,distance_to_J1_um=distj1,distance_to_J2_um=distj2,
            last_snapshot_x_um=position[0]*1e6,last_snapshot_y_um=position[1]*1e6,last_snapshot_z_um=position[2]*1e6,
            snapshot_to_final_distance_m=float(np.linalg.norm(position-tail[1:4])),
            local_fluid_speed_mm_s=float(np.linalg.norm(sample.velocity_m_s)*1e3),
            unconstrained_hydro_speed_mm_s=float(np.linalg.norm(solution.unconstrained[:3])*1e3),
            constrained_particle_speed_mm_s=float(np.linalg.norm(solution.velocity[:3])*1e3),
            final_angular_speed_s_inv=float(np.linalg.norm(solution.velocity[3:])),
            pressure_Pa=float(sample.pressure_pa),tetra_id=int(sample.tetra_id),
            final_contacts=len(current_contacts),independent_constraint_rank=solution.record['contact_redundancy']['rank_after'],
            minimum_contact_multiplier=float(np.min(solution.record['multipliers'])),
            max_contact_normal_angle_deg=angle,geometry_gap_nm=plane.gap_m*1e9,handoff_lower_nm=lower*1e9,
            last64_position_span_m=span,stationary_position_budget_m=16*wall.roundoff_m,
            trailing_identical_nominal_steps=trailing,near_zero_translation_onset_ms=onset*1e3,
            diagnostic_stationary_tail_ms=(tail[0]-onset)*1e3,downstream_direction_test=lp['status'],
            max_downstream_projection=lp['maximum_forward_projection'],no_downstream_certificate_residual=lp['dual_certificate_residual'],
            nonzero_escape_direction_exists=lp['any_nonzero_feasible_direction'],
            static_translation_replay_error_m_s=verr,static_hydro_replay_error_m_s=herr,
            conditioning=solution.record['condition_estimate'],backward_residual=solution.record['relative_backward_residual'])
        stops.append(dict(base,**extra))
        static.append(dict(particle_id=pid,translation_replay_error_m_s=verr,hydro_replay_error_m_s=herr,
            identical_contact_ids=True,all_32_recorded_solves_pass=True,original_audit_window_matches=True,
            maximum_forward_projection=lp['maximum_forward_projection'],lp_status=lp['status'],
            dual_certificate_residual=lp['dual_certificate_residual'],escape_direction_exists=lp['any_nonzero_feasible_direction']))
        if len(stops)%40==0:print('stops analyzed',len(stops),flush=True)
    assert len(stops)==185 and not any(safety.values())
    csvout('all_tracks.csv',allrows);csvout('stopped_tracks.csv',stops);csvout('contact_geometry.csv',contacts_csv);csvout('terminal_static_checks.csv',static)
    regionrows=[]
    for region,count in collections.Counter(r['region'] for r in stops).items():
        rr=[r for r in stops if r['region']==region]
        regionrows.append(dict(region=region,count=count,diameter_min_um=min(r['diameter_um'] for r in rr),
            diameter_max_um=max(r['diameter_um'] for r in rr),termination_min_ms=min(r['termination_age_ms'] for r in rr),
            termination_max_ms=max(r['termination_age_ms'] for r in rr)))
    csvout('regions.csv',regionrows)
    # Stable examples: smallest, median and largest stopped bubble; keep ID 2 for an easily inspected original record.
    ds=sorted(stops,key=lambda r:r['diameter_um']);example_ids=list(dict.fromkeys([ds[0]['particle_id'],ds[len(ds)//2]['particle_id'],ds[-1]['particle_id'],2]))
    example_rows=[]
    for pid in example_ids:
        s=stopped_raw[pid];us=field.sample_many(s[:,1:4]);assert us.inside_lumen.all()
        for i,v in enumerate(s):example_rows.append(dict(particle_id=pid,age_ms=v[0]*1e3,x_um=v[1]*1e6,y_um=v[2]*1e6,z_um=v[3]*1e6,
            particle_speed_mm_s=float(np.linalg.norm(v[4:7])*1e3),fluid_speed_mm_s=float(np.linalg.norm(us.velocity_m_s[i])*1e3),
            wall_gap_nm=v[14]*1e9,accepted_dt_ms=v[18]*1e3))
    csvout('example_timeseries.csv',example_rows)
    # Emulate the existing four movie views' shared endpoint policy, without rendering or altering them.
    movie=read(RUN/'data/render_manifest.json');age=float(movie['age_end_s'])
    held=0
    for r in stops:
        s=stopped_raw[r['particle_id']]
        p=np.array([np.interp(age,s[:,0],s[:,k]) for k in [1,2,3]])
        assert np.array_equal(p,s[-1,1:4]);held+=1
    shape_stats={}
    for name in ['COMPLETED','SUPPORTED_STATIONARY']:
        rr=[r for r in allrows if r['status']==name]
        shape_stats[name]={k:distribution([r[k] for r in rr]) for k in ['diameter_um','termination_age_ms','path_length_um']}
    output=dict(all_pass=True,flow_sha256=summary['flow_sha256'],dt_s=.0005,total=1500,status_counts=dict(collections.Counter(r['status'] for r in rows)),
        outlet_counts=dict(collections.Counter(r['outlet'] for r in rows if r['outlet'])),
        stop_reason_counts=dict(collections.Counter(r['failure_detail'] for r in stops)),stopped_fraction=len(stops)/1500,
        terminal_regions=regionrows,source_case='ROI-only-balanced-pressure-v1',protected_sources_verified=len(sources),track_files_hash_verified=verified_files,
        source_data_identity_verified=True,wall_array_identity_verified=True,original_audit_windows_rechecked=same_last_windows,
        groups=shape_stats,stop_statistics={k:distribution([r[k] for r in stops]) for k in ['local_fluid_speed_mm_s','unconstrained_hydro_speed_mm_s',
            'constrained_particle_speed_mm_s','final_angular_speed_s_inv','geometry_gap_nm','final_contacts','termination_age_ms',
            'static_translation_replay_error_m_s','static_hydro_replay_error_m_s','max_contact_normal_angle_deg','conditioning','backward_residual']},
        forward_feasibility_counts=dict(collections.Counter(r['downstream_direction_test'] for r in stops)),
        nonzero_escape_direction_exists_count=sum(r['nonzero_escape_direction_exists'] for r in stops),
        stop_center_bounds_um=[np.min([r['final_position_m'] for r in rows if r['status']=='SUPPORTED_STATIONARY'],axis=0)*1e6,
            np.max([r['final_position_m'] for r in rows if r['status']=='SUPPORTED_STATIONARY'],axis=0)*1e6],
        wall_roundoff_m=wall.roundoff_m,stagnation_budget_m=16*wall.roundoff_m,example_ids=example_ids,
        movie_policy=dict(final_age_s=age,stopped_endpoints_held=held,original_video_seconds=movie['video_seconds'],
            source='render_results.py:172-173; additional_views/render_views.py:60-65',stationary_endpoint_extension_is_visual_only=True),
        safety_totals=dict(safety),new_CFD_runs=0,new_trajectory_integrations=0,terminal_static_solves=185,
        geometric_LP_scope='Instantaneous fixed-size sphere constraints, not global passage feasibility or physiological trapping proof',
        elapsed_s=time.time()-start)
    dump('data/analysis_summary.json',output);dump('data/source_hashes.json',LOCK)
    for p,h in LOCK.items():assert sha(p)==h,p
    print(json.dumps(output,ensure_ascii=False,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()),flush=True)

if __name__=='__main__':main()
