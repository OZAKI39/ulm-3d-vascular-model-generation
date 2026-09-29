from pathlib import Path
import argparse,json,sys,time
import numpy as np
import torch

ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.population_v1.provenance import write_json
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_gap import wall_gap
from particle_3d.nearfield_handoff import handoff_constraints,first_handoff_event,wall_handoff_certificate,partitioned_wall_handoff_certificate
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.open_boundary_rim import rim_witness

root=Path(__file__).resolve().parents[1];env=load_environment(json.loads((root/'data/environment.json').read_text()))
print('ENVIRONMENT_LOADED',flush=True);gpu=GPUWall(env.wall);policy=NearFieldRegularizationV1()
rng=np.random.default_rng(2026092907);ids=rng.integers(0,len(env.wall.triangles),800)
radius=rng.uniform(.375e-6,2e-6,len(ids));gap=np.tile([2e-9,10e-9,1e-6,50e-6],len(ids)//4)
x=env.wall.triangle_centers[ids]+(radius+gap)[:,None]*env.wall.normal_in[ids]
# Add edge/vertex configurations, including actual open rim topology.
for feature in ['edge','vertex']:
    points=env.wall.triangles[ids[:100],:2].mean(1) if feature=='edge' else env.wall.triangles[ids[:100],0]
    x=np.vstack([x,points+(radius[:100]+gap[:100])[:,None]*env.wall.normal_in[ids[:100]]])
radii=np.r_[radius,radius[:100],radius[:100]]
tensor=torch.tensor(x,device='cuda');rt=torch.tensor(radii,device='cuda')
ni,value,contacts,normal,count=gpu.query(tensor,rt);ni=ni.cpu().numpy();value=value.cpu().numpy();contacts=contacts.cpu().numpy();normal=normal.cpu().numpy();count=count.cpu().numpy()
errors=dict(gap_m=0.,wall_point_m=0.,normal=0.);mismatch=[];contact_mismatch=[];rim_mismatch=[];normal_bad=[]
start=time.perf_counter()
for k,(point,r) in enumerate(zip(x,radii)):
    s=Sphere(point,r);ref=wall_gap(s,env.wall)
    if ni[k]!=ref.wall_triangle_id:mismatch.append(dict(index=k,cpu=ref.wall_triangle_id,gpu=int(ni[k]),gap_difference_m=float(value[k,0]-ref.gap_m)))
    errors['gap_m']=max(errors['gap_m'],abs(value[k,0]-ref.gap_m))
    errors['wall_point_m']=max(errors['wall_point_m'],float(np.max(abs(value[k,1:4]-ref.wall_point_m))))
    ne=float(np.max(abs(value[k,7:10]-ref.normal_inward)));errors['normal']=max(errors['normal'],ne)
    if ne>1e-9:normal_bad.append(k)
    rim=bool(rim_witness(s,env.wall,ref));
    if rim!=bool(value[k,13]):rim_mismatch.append(k)
    if k<100:
        cs,_=handoff_constraints({1:s},env.wall,policy);ref_ids=sorted(c.canonical_id[3] for c in cs if c.canonical_id[0]=='CONTINUUM_HANDOFF_CONTACT')
        got=contacts[k,:count[k]].tolist()
        if got!=ref_ids:contact_mismatch.append(dict(index=k,cpu=ref_ids,gpu=got))
print('NEAREST_AND_CONTACT_PARITY',time.perf_counter()-start,flush=True)
# A deterministic selection of paths that start above the inherited handoff.
select=[k for k in range(len(x)) if value[k,0]>=2e-9-env.wall.roundoff_m][:60]
px=x[select];pr=radii[select];vel=rng.normal(0,.03,(len(select),3));dt=.0005
path=gpu.path(torch.tensor(px,device='cuda'),torch.tensor(vel,device='cuda'),torch.tensor(pr,device='cuda'),dt).cpu().numpy()
root_errors=[];false_certificates=[];missed_roots=[]
for k,(p,r,v) in enumerate(zip(px,pr,vel)):
    shape=Sphere(p,r);hit=first_handoff_event({1:shape},{1:v},dt,env.wall,policy)
    if hit is not None:
        if path[k,1]<0:missed_roots.append(k)
        else:root_errors.append(abs(hit['time_from_start_s']-path[k,0]))
    if path[k,2]==1:
        safe,proof=wall_handoff_certificate(shape,shape.moved(p+dt*v),env.wall,2e-9)
        if not safe:safe,proof=partitioned_wall_handoff_certificate(shape,shape.moved(p+dt*v),env.wall,2e-9)
        if not safe:false_certificates.append(dict(index=k,cpu=proof))
record=dict(PASS=not mismatch and not contact_mismatch and not rim_mismatch and not normal_bad and not missed_roots and not false_certificates
    and errors['gap_m']<=8*env.wall.roundoff_m and max(root_errors,default=0.)<1e-10,
    phase=7,scope='GPU_WALL_GEOMETRY_AND_CONTINUOUS_PATH_ONLY_NOT_FULL_DYNAMICS',flow_sha256=env.flow_sha256,
    query_count=len(x),errors=errors,nearest_triangle_mismatches=mismatch[:20],contact_mismatches=contact_mismatch,
    rim_mismatches=rim_mismatch,normal_mismatches=normal_bad,wall_roundoff_m=env.wall.roundoff_m,
    path_count=len(select),maximum_root_time_error_s=max(root_errors,default=0.),missed_roots=missed_roots,
    false_positive_certificates=false_certificates,certified_on_gpu=int(sum(path[:,2]==1)),
    needs_finite_triangle_partition_or_subdivision=int(sum(path[:,2]!=1)))
write_json(root/'validation/phase7_real_wall_gpu.json',record);print((root/'validation/phase7_real_wall_gpu.json').read_text(),flush=True)
if not record['PASS']:raise SystemExit(1)
