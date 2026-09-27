"""Does a C representation skip a WALL conflict? Audit, never change births.

The swept volume of a fixed sphere on each saved linear path segment is a
capsule. Existing exact finite-triangle capsule gaps certify the whole segment.
Open caps are excluded from solid-wall constraints throughout this transition.
"""
from pathlib import Path
import argparse,gzip,json,multiprocessing as mp,socket,time
import numpy as np
from .particle82a_admission import context
from .particle82a_geometry import lower_gap
from .particle_shapes import Capsule,Sphere
from .wall_gap import wall_gap
from .particle82_provenance import atomic_json,require_remote


def swept_entry(points,radius,wall):
    lower=float(lower_gap(radius));ro=wall.roundoff_m
    initial=wall_gap(Sphere(points[0],radius),wall).gap_m
    if initial< -ro:
        return dict(wall_geometrically_clear=False,handoff_clear=False,
            first_failure='WALL_INTERSECTION_AT_ANCHOR',
            anchor_wall_gap_m=initial,minimum_evaluated_gap_m=initial,
            full_sweep_evaluated=False,segment_checks=0,roundoff_m=ro)
    smallest=initial;checks=0
    for a,b in zip(points[:-1],points[1:]):
        delta=b-a;length=float(np.linalg.norm(delta))
        if length==0:continue
        swept=Capsule((a+b)/2,delta/length,radius,length)
        gap=wall_gap(swept,wall).gap_m;smallest=min(smallest,gap);checks+=1
    return dict(wall_geometrically_clear=bool(smallest>=-ro),handoff_clear=bool(smallest>=lower-ro),
        first_failure=None if smallest>=lower-ro else 'WALL_INTERSECTION_ON_ENTRY_PATH' if smallest< -ro else 'NEARFIELD_HANDOFF_VIOLATION_ON_ENTRY_PATH',
        anchor_wall_gap_m=initial,minimum_evaluated_gap_m=smallest,full_sweep_evaluated=True,
        segment_checks=checks,roundoff_m=ro)


def job(file):
    file=Path(file);c=context()
    with gzip.open(file,'rt') as f:rows=json.load(f)['rows']
    raw=np.load(file.with_name(file.name.replace('.json.gz','.npz')));result=[]
    for i,row in enumerate(rows):
        event=row['event'];method=row['methods']['C']
        if not method['accepted']:continue
        path=raw['paths'][raw['offsets'][i]:raw['offsets'][i+1],1:]
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
        keep=path[arc<method['s_birth_m']]
        points=np.vstack([keep,np.asarray(method['birth_center_m'])])
        audit=swept_entry(points,method['radius_m'],c.env.wall)
        result.append(dict(event_id=event['event_id'],anchor_basin=event['point_tracer_basin'],
            birth_basin=method['birth_point_basin'],s_birth_m=method['s_birth_m'],radius_m=method['radius_m'],**audit))
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--admission',required=True);p.add_argument('--report',required=True)
    p.add_argument('--provenance',required=True);p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    prov=json.loads(Path(a.provenance).read_text());require_remote(prov,prov['hostname']);context();start=time.time()
    files=sorted((Path(a.admission)/'events').glob('events_*.json.gz'));rows=[]
    with mp.get_context('fork').Pool(a.workers) as pool:
        for result in pool.imap_unordered(job,map(str,files),chunksize=1):rows.extend(result)
    rows.sort(key=lambda r:r['event_id']);summary={}
    from collections import Counter
    for b in ['ALL','OUTLET_01','OUTLET_02','OUTLET_03','UNRESOLVED_POINT_PATH']:
        selected=rows if b=='ALL' else [r for r in rows if r['anchor_basin']==b]
        summary[b]=dict(C_births_found=len(selected),wall_geometrically_clear=sum(r['wall_geometrically_clear'] for r in selected),
            full_sweep_handoff_certified=sum(r['handoff_clear'] for r in selected),
            failure_counts=dict(Counter(r['first_failure'] for r in selected if r['first_failure'])))
    atomic_json(Path(a.report)/'METHOD_C_ENTRY_SWEEP_AUDIT.json',dict(summary=summary,rows=rows,
        birth_results_unchanged=True,geometry_only_not_a_force_or_entry_dynamics_model=True,
        certificate='EXISTING_EXACT_CAPSULE_VS_ORIGINAL_WALL_TRIANGLES_FOR_EVERY_SAVED_LINEAR_SEGMENT',
        open_caps_are_not_solid=True,REMOTE_SERVER_COMPUTE=True,hostname=socket.gethostname(),
        source_commit=prov['source_git_commit'],seconds=time.time()-start))


if __name__=='__main__':main()
