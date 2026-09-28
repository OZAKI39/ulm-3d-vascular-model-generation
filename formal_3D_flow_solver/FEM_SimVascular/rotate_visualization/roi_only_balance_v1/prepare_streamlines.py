"""Finite-element streamlines of the new 2 mm/s field, with certified cap exits.

Reuse the existing double-precision P1 RK23 point tracer. No particle forces,
radius, seed offsets, target-outlet steering, or post-hoc path extension.
"""
from pathlib import Path
from types import SimpleNamespace
from collections import Counter
import argparse,csv,hashlib,json,socket,sys,time
import numpy as np
import pyvista as pv

ROOT=Path(__file__).resolve().parent
CASE=ROOT/'input_data';OUT=CASE/'streamlines'
PARTICLE_SOURCE=Path('/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/src')
QUOTAS={'OUTLET_01':16,'OUTLET_02':56,'OUTLET_03':24}
STEP=.15e-6;ERROR=1e-11;SEED=2026092203


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def setup(source=PARTICLE_SOURCE):
    sys.path.insert(0,str(source))
    from particle_3d.field import FrozenFEMField
    from particle_3d.inlet_flux import InletFluxSampler
    from particle_3d.validation_boundary import ValidationBoundaryClassifier
    from particle_3d.particle82_point_native import NativePointTracer
    a=np.load(CASE/'frozen_flow/flow_arrays_si.npz')
    field=FrozenFEMField(a['points_m'],a['tetra'],a['velocity_m_s'],a['pressure_pa'])
    boundaries={role:pv.read(CASE/'solver_mesh/mesh-surfaces'/f'{role}.vtp') for role in ['INLET',*QUOTAS]}
    classifier=ValidationBoundaryClassifier({k:v for k,v in boundaries.items() if k in QUOTAS})
    env=SimpleNamespace(field=field,boundaries=boundaries,classifier=classifier)
    native=NativePointTracer(env)
    inlet=boundaries['INLET'];ids=np.asarray(inlet['GlobalNodeID'],int)[inlet.faces.reshape(-1,4)[:,1:]]-1
    triangles=a['points_m'][ids];normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    velocity=a['velocity_m_s'][ids]
    normals[np.einsum('ni,ni->n',normals,velocity.mean(1))<0]*=-1
    flux=np.einsum('nki,ni->nk',velocity,normals)
    sampler=InletFluxSampler(triangles,flux)
    measured=json.loads((CASE/'reports/physics_validation_H0.json').read_text())['measurements']['Q_in_m3_s']
    assert np.isclose(sampler.Q_m3_s,measured,rtol=1e-12,atol=0)
    return env,native,sampler


def spatial_subset(records,quota):
    """Farthest-point display selection spreads actual endpoints over each cap."""
    if len(records)<quota:raise ValueError('Insufficient naturally completed streamlines')
    points=np.array([r['path'][-1,1:] for r in records]);center=points.mean(0)
    indices=[int(np.argmin(np.linalg.norm(points-center,axis=1)))]
    nearest=np.linalg.norm(points-points[indices[0]],axis=1)
    while len(indices)<quota:
        nearest[indices]=-1;index=int(np.argmax(nearest));indices.append(index)
        nearest=np.minimum(nearest,np.linalg.norm(points-points[index],axis=1))
    return [records[k] for k in indices]


def resampled_path(path,count=1001):
    s=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(path[:,1:],axis=0),axis=1))]
    return np.column_stack([np.interp(np.linspace(0,s[-1],count),s,path[:,k]) for k in [1,2,3]])


def main():
    p=argparse.ArgumentParser();p.add_argument('--particle-source',type=Path,default=PARTICLE_SOURCE)
    p.add_argument('--candidates',type=int,default=768);a=p.parse_args()
    OUT.mkdir(exist_ok=True);(OUT/'data').mkdir(exist_ok=True)
    if (OUT/'COMPUTE_VALIDATION.json').exists():raise RuntimeError('Existing completed streamlines: render/review saved data; do not recompute')
    paths=list((CASE/'frozen_flow').rglob('*'))+list((CASE/'solver_mesh').rglob('*'))
    paths+=list((CASE/'figures').glob('*'))+list((CASE/'animations').glob('*'))
    paths+=[CASE/'reports/physics_validation_H0.json',CASE/'reports/physics_validation_H0.json',CASE/'reports/render_manifest.json']
    lock={str(p.relative_to(CASE)):sha(p) for p in paths if p.is_file()}
    dump(OUT/'SOURCE_LOCK.json',lock);start=time.time();env,native,sampler=setup(a.particle_source)
    seeds,triangles=sampler.sample(np.random.default_rng(SEED),a.candidates)
    records=[];counts=Counter();candidates={}
    for i,position in enumerate(seeds):
        trace=native.trace(position,step_m=STEP,error=ERROR,horizon_m=500e-6)
        row=dict(seed_id=i,inlet_triangle=int(triangles[i]),outlet=trace['outlet'],end_reason=trace['end_reason'],path=trace['path'])
        records.append(row);counts[row['end_reason']]+=1;candidates[f'seed_{i:04d}']=trace['path']
        if (i+1)%96==0:print('STREAMLINE_CANDIDATES',i+1,'/',a.candidates,dict(counts),flush=True)
    np.savez_compressed(OUT/'data/candidates.npz',**candidates)
    dump(OUT/'data/candidate_ledger.json',[{k:v for k,v in r.items() if k!='path'} for r in records])
    selected=[]
    for role,quota in QUOTAS.items():selected.extend(spatial_subset([r for r in records if r['outlet']==role],quota))
    point_arrays=[];velocity_arrays=[];cells=[];offset=0;catalog=[];archive={};refined={}
    for index,r in enumerate(selected):
        path=r['path'];positions=np.ascontiguousarray(path[:,1:]);sampled=np.empty_like(positions)
        cell=env.field.locate(positions[0])[0]
        for k,position in enumerate(positions):
            vector=np.empty(3);cell=native.lib.p82_sample(position.ctypes.data,cell,vector.ctypes.data)
            assert cell>=0,'Saved streamline point outside original P1 mesh'
            sampled[k]=vector
        # This is a new integration of the same seed, not a resampled display curve.
        fine=native.trace(positions[0],step_m=STEP/2,error=ERROR/4,horizon_m=500e-6)
        assert fine['outlet']==r['outlet'],'Outlet changed on spatial-step refinement'
        error=float(np.linalg.norm(resampled_path(path)-resampled_path(fine['path']),axis=1).max())
        endpoint_error=float(np.linalg.norm(path[-1,1:]-fine['path'][-1,1:]))
        assert error<.05e-6 and endpoint_error<.05e-6,'Display streamline has not converged geometrically'
        hit=env.classifier.first_event(path[-2,1:],path[-1,1:]);assert hit and hit.role==r['outlet']
        assert np.isfinite(path).all() and np.all(np.diff(path[:,0])>0)
        point_arrays.append(positions);velocity_arrays.append(sampled)
        cells.extend([len(path),*range(offset,offset+len(path))]);offset+=len(path)
        archive[f'line_{index:03d}']=path;archive[f'velocity_{index:03d}']=sampled;refined[f'line_{index:03d}']=fine['path']
        catalog.append(dict(line_id=index,seed_id=r['seed_id'],outlet=r['outlet'],inlet_triangle=r['inlet_triangle'],
            endpoint_triangle=hit.role_triangle_id,points=len(path),transit_time_s=float(path[-1,0]),
            length_um=float(np.linalg.norm(np.diff(positions,axis=0),axis=1).sum()*1e6),
            path_refinement_error_um=error*1e6,endpoint_refinement_error_um=endpoint_error*1e6,
            initial_position_m=positions[0].tolist(),terminal_position_m=positions[-1].tolist()))
        if (index+1)%24==0:print('STREAMLINE_VALIDATED',index+1,'/',len(selected),flush=True)
    mesh=pv.PolyData(np.vstack(point_arrays),lines=np.array(cells));mesh['Velocity_m_s']=np.vstack(velocity_arrays)
    mesh['Speed_mm_s']=np.linalg.norm(mesh['Velocity_m_s'],axis=1)*1000
    mesh.cell_data['Outlet_id']=np.array([int(r['outlet'][-2:]) for r in selected]);mesh.cell_data['Line_id']=np.arange(len(selected))
    mesh.save(OUT/'data/streamlines_si.vtp');np.savez_compressed(OUT/'data/selected_paths.npz',**archive)
    np.savez_compressed(OUT/'data/refined_paths.npz',**refined);dump(OUT/'data/catalog.json',catalog)
    with (OUT/'data/catalog.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(catalog[0]));writer.writeheader();writer.writerows(catalog)
    changed=[p for p,h in lock.items() if sha(CASE/p)!=h];assert not changed
    dependencies={str(p.relative_to(a.particle_source)):sha(p) for p in a.particle_source.rglob('*') if p.suffix in ['.py','.cpp']}
    result=dict(all_pass=True,field_case='mean-2p0-mmps-A-ROI-only-balanced-pressure-v1',inlet_mean_mm_s=2.,
        source_field_sha256=sha(CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu'),
        source_npz_sha256=sha(CASE/'frozen_flow/flow_arrays_si.npz'),
        seed=SEED,candidate_count=len(records),candidate_outcomes=dict(counts),selected_count=len(selected),
        selected_outlet_counts=QUOTAS,selection='OUTLET_STRATIFIED_FARTHEST_ENDPOINT_SPACING_FOR_VISUALIZATION',
        line_counts_are_not_flow_fractions=True,all_seeds_on_original_inlet=True,
        all_selected_endpoints_reclassified_on_actual_outlet_triangles=True,
        all_selected_samples_inside_original_tetra_mesh=True,all_selected_refined_outlet_labels_match=True,
        maximum_refined_path_difference_um=max(r['path_refinement_error_um'] for r in catalog),
        maximum_refined_endpoint_difference_um=max(r['endpoint_refinement_error_um'] for r in catalog),
        method='DOUBLE_PRECISION_P1_RK23_NATIVE',maximum_spatial_step_m=STEP,local_position_error_tolerance_m=ERROR,
        terminal_step_rule='ACTUAL_OPEN_CAP_INTERSECTION_WITH_LOCAL_EULER_STEP_WITHIN_2_PERCENT_SPATIAL_STEP',
        terminal_approximation_maximum_length_nm=STEP*.02*1e9,
        no_target_outlet_steering=True,no_artificial_path_extension=True,no_spline_smoothing=True,
        physical_role='STEADY_FLOW_POINT_STREAMLINES_NOT_FINITE_SIZE_MICROBUBBLE_TRAJECTORIES',
        original_inputs_and_visuals_unchanged=True,unchanged_file_count=len(lock),
        dependency_root=str(a.particle_source),dependency_sha256=dependencies,
        outputs_sha256={str(p.relative_to(OUT)):sha(p) for p in (OUT/'data').iterdir() if p.is_file()},
        hostname=socket.gethostname(),elapsed_seconds=time.time()-start,script_sha256=sha(__file__))
    dump(OUT/'COMPUTE_VALIDATION.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['dependency_sha256','outputs_sha256']},indent=2),flush=True)


if __name__=='__main__':main()
