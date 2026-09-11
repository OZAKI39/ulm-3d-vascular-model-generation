"""CPU inspection of the authorized short controls; never launches a solver."""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from py_scripts.fluid_physics.common import now, sha256_file, write_json
from py_scripts.single_rbc_benchmark.analysis import read_csv
from py_scripts.single_rbc_benchmark.physics import read_off, geometry
from py_scripts.single_rbc_benchmark.quality import mesh_intersections
from .geometry_checks import checked_mesh, membership, independent_intersections, require_same_frame


def inspect_run(directory, output):
    start=time.perf_counter();d=Path(directory).resolve();out=Path(output).resolve()
    out.mkdir(parents=True,exist_ok=True)
    result=dict(task=d.name,source_directory=str(d),recorded_at=now(),scope='Saved short-control states only; not a full-range fix or material acceptance',qualified_speedup=None)
    execution=json.loads((d/'execution.json').read_text());result['execution']=execution
    completion=d/'completion.json';result['completion']=json.loads(completion.read_text()) if completion.exists() else None
    spec=json.loads((d/'spec.json').read_text());c=spec['config'];reference,faces=read_off(spec['mesh'])
    raw=read_csv(d/'vertices.csv') if (d/'vertices.csv').exists() else None
    if raw is None or not len(raw):
        result.update(status='NO_SAVED_MEMBRANE_STATE',cpu_analysis_s=time.perf_counter()-start)
        write_json(out/'review.json',result);return result
    moments=read_csv(d/'moments.csv');ref=geometry(reference,faces)
    edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
    lengths0=np.linalg.norm(reference[edges[:,1]]-reference[edges[:,0]],axis=1)
    first_member=first_cross=None;rows=[];input_files={d/'spec.json',d/'execution.json',d/'vertices.csv',d/'moments.csv',Path(spec['mesh'])}
    keys=list(dict.fromkeys(zip(raw['phase'],raw['step'])))
    total_probes=total_mismatch=total_uncertain=total_ambiguous=0
    native_ids_verified=True
    for frame_index,(phase,step) in enumerate(keys):
        frame=raw[(raw['phase']==phase)&(raw['step']==step)]
        v,g=checked_mesh(frame['vertex'],np.column_stack([frame[k] for k in ('x','y','z')]),faces,c['geometry']['periodic_length'])
        identity=dict(phase=str(phase),step=int(step),time_star=float(frame[0]['time_star']),strain=float(frame[0]['strain']))
        assert len(np.unique(frame['strain']))==len(np.unique(frame['time_star']))==1
        done=int(step)-(spec['prep_steps'] if phase=='shear' else 0)
        p=d/f'fluid_probes_{phase}_{done:08d}.npz';state=d/f'membrane_state_{phase}_{done:08d}.npz';input_files.add(p)
        with np.load(p,allow_pickle=False) as z:
            require_same_frame(phase,step,str(z['phase'].item()),int(z['global_step']))
            assert abs(float(z['strain'])-identity['strain'])<1e-12
            cloud=z['probes'].copy();ids=z['membrane_vertex_ids'].copy()
        assert np.array_equal(ids,np.arange(len(reference))), 'NATIVE_MEMBRANE_ID_CONNECTIVITY_UNVERIFIED'
        velocity_scale=None
        if state.exists():
            input_files.add(state)
            with np.load(state,allow_pickle=False) as z:
                require_same_frame(phase,step,str(z['phase'].item()),int(z['global_step']))
                assert np.array_equal(z['ids'],ids)
                assert np.allclose(z['positions'],np.column_stack([frame[k] for k in ('x','y','z')]),rtol=0,atol=1e-9)
                velocity_scale=float(np.linalg.norm(z['velocities'],axis=1).max()*spec['dt']/lengths0.min())
        else:native_ids_verified=False
        points=cloud[:,2:5].copy();L=c['geometry']['periodic_length']
        points[:,:2]+=L*np.rint((v.mean(axis=0)[:2]-points[:,:2])/L)
        check=membership(points,cloud[:,0]==1,v,faces,c['repair']['containment_band'])
        n=int(check['confirmed_mismatch'].sum());uncertain=int(check['uncertain'].sum());ambiguous=int(check['ambiguous'].sum())
        total_probes+=len(points);total_mismatch+=n;total_uncertain+=uncertain;total_ambiguous+=ambiguous
        if n and first_member is None:
            i=int(np.flatnonzero(check['confirmed_mismatch'])[0]);tri=int(check['nearest_triangle'][i])
            first_member=dict(**identity,particle_id=int(cloud[i,1]),membership=int(cloud[i,0]),position=points[i].tolist(),distance_to_membrane=float(check['distance'][i]),triangle=tri,triangle_vertex_ids=faces[tri].tolist(),triangle_coordinates=v[faces[tri]].tolist(),confirmed_in_frame=n)
        legacy=mesh_intersections(v,faces);independent=None
        if legacy['nonadjacent_intersections'] or frame_index==len(keys)-1:
            independent=independent_intersections(v,faces,c['repair']['intersection_tolerance'])
            if independent['confirmed_count'] and first_cross is None:first_cross=dict(**identity,**independent)
        m=moments[(moments['phase']==phase)&(moments['step']==step)];assert len(m)==1
        lengths=np.linalg.norm(v[edges[:,1]]-v[edges[:,0]],axis=1)
        rows.append(dict(**identity,area=g['area'],volume=g['volume'],D=g['D'],theta_deg=g['theta_deg'],area_relative_drift=abs(g['area']/ref['area']-1),volume_relative_drift=abs(g['volume']/ref['volume']-1),minimum_face_area=g['minimum_face_area'],max_wlc_extension_fraction=float(np.max(lengths/lengths0)*c['mirheo_membrane']['x0']),membrane_speed_dt_over_reference_min_edge=velocity_scale,confirmed_member_point_frames=n,near_surface_uncertain=uncertain,ambiguous_rays=ambiguous,legacy_intersection_count=legacy['nonadjacent_intersections'],independent_intersection_count=None if independent is None else independent['confirmed_count'],fluid_N=int(m[0]['N']),inner_N=int(m[0]['inner_N']),temperature=float(m[0]['temperature']),fluid_max_speed=float(m[0]['max_speed']),wall_crossings=int(m[0]['wall_crossings'])))
    with (out/'frame_metrics.csv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    actual=result['completion'];finished=bool(actual and actual.get('completed') and actual.get('actual_steps')==spec['prep_steps']+spec['steps'])
    screening=dict(same_phase_native_ids_and_positions=native_ids_verified,no_confirmed_probe_discrepancy=total_mismatch==0,no_confirmed_self_intersection=first_cross is None,area_within_original_2pct=max(x['area_relative_drift'] for x in rows)<=c['screen']['area_relative_drift'],volume_within_original_2pct=max(x['volume_relative_drift'] for x in rows)<=c['screen']['volume_relative_drift'],constant_fluid_N=len({x['fluid_N'] for x in rows})==1,no_saved_wall_crossings=all(x['wall_crossings']==0 for x in rows),completed_planned_short_protocol=finished)
    logs={}
    for p in d.glob('*00000.log'):
        input_files.add(p)
        failures=[x for x in p.read_text(errors='replace').splitlines() if 'too many triangle collision candidates' in x]
        if failures:logs[p.name]=failures
    result.update(status='SAVED_SHORT_CHECKS_PASS' if all(screening.values()) else 'SHORT_CONTROL_HAS_UNRESOLVED_QUALITY_FAILURE',short_screen=screening,frames=len(rows),last_saved_frame=rows[-1],first_confirmed_membership_discrepancy=first_member,first_confirmed_self_intersection=first_cross,probe_summary=dict(tested_point_frames=total_probes,confirmed_mismatched_point_frames=total_mismatch,near_surface_uncertain_point_frames=total_uncertain,ambiguous_ray_point_frames=total_ambiguous,scope='Fixed identity sample, not whole-fluid leakage or continuous-time proof'),native_error_lines=logs,limitations=['A short pass does not cover Gamma=4 or the original shear-failure range.','Velocity times dt is a diagnostic scale, not measured one-step oldPositions displacement.','No material, force-feedback or general impermeability qualification follows from this saved-state screen.'],source_sha256={str(p):sha256_file(p) for p in sorted(input_files)},cpu_analysis_s=time.perf_counter()-start)
    write_json(out/'review.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    r=inspect_run(a.run,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('source_sha256','execution')},ensure_ascii=False,allow_nan=False,indent=2))


if __name__=='__main__':main()
