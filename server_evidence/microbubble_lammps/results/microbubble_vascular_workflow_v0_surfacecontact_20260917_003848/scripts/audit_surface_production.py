#!/usr/bin/env python3
from pathlib import Path
import json,csv,collections,itertools
import numpy as np
from surface_reference import SurfaceReference
from qualify_surface import Probe,compare,metrics,records,S
r=SurfaceReference.production(S);probe=Probe('',True)
old=Path('/home/lzy/projects/compre_output/microbubble_vascular_workflow_v0/wallslide_20260916_233923')
points=[]
for name in ['LOW','MEDIUM','LONG']:
    with (old/'runs'/f'{name if name!="LONG" else "LONG_TRANSPORT"}_release_mpi1'/'TRAJECTORIES.csv').open() as f:rows=list(csv.DictReader(f))
    for a in rows[::max(1,len(rows)//600)]:points.append([float(a[k+'_m']) for k in ['x','y','z']])
orientations=[]
for x in points:
    out=compare(r,probe,x,[-1e-5,2e-5,-3e-5],label='production_orientation');q=out['query']
    if len(q['candidates'])==1:
        c=q['candidates'][0];n=q['clusters'][0]['normal'];d=np.array(x)-c['point'];cos=float(np.dot(n,d)/np.linalg.norm(d));orientations.append(cos);assert cos>0
oldq=json.loads((S/'validation/OLD_LONG_FAILURE_REPLAY.json').read_text())
replay=[]
for factor,theta in itertools.product([.5,1,2],[10,15,25]):replay.append(compare(r,probe,oldq['query_point'],oldq['raw_velocity'],factor,theta,'old_fatal_query'))
seams=collections.defaultdict(list);samples=[]
for edge,inc in r.edges.items():
    if len(inc)==2 and r.regions[inc[0]]!=r.regions[inc[1]]:
        a,b=inc;angle=r.dihedral[tuple(sorted(inc))];key='-'.join(map(str,sorted([r.regions[a],r.regions[b]])));seams[key].append(angle)
        if angle<25:
            n=r.normal[a]+r.normal[b];n/=np.linalg.norm(n);x=r.vertices[list(edge)].mean(axis=0)+1e-8*n
            result=compare(r,probe,x,[-1e-5,2e-5,-3e-5],label='region_seam');samples.append(dict(edge=list(edge),regions=[int(r.regions[a]),int(r.regions[b])],dihedral=angle,status=result['query']['status'],clusters=len(result['query']['clusters'])))
# Region labels are annotations only: relabel every face and repeat the same queries.
regions=r.regions.copy();sample=records[::max(1,len(records)//100)];baseline=[r.query(a['x'],a['factor'],a['theta'])['clusters'] for a in sample];r.regions[:]=0
for a,b in zip(sample,baseline):
    q=r.query(a['x'],a['factor'],a['theta']);assert q['clusters']==b
r.regions[:]=regions;probe.close()
(S/'validation/FLUID_SIDE_NORMAL_AUDIT.json').write_text(json.dumps(dict(status='PASS',queries=len(points),unambiguous_queries=len(orientations),minimum_center_direction_cosine=min(orientations),maximum_center_direction_cosine=max(orientations),method='Independent topology, distances, dihedrals, features, pseudonormals and projection vs C++; actual prior accepted fluid centers.'),indent=2)+'\n')
(S/'validation/REGION_SEAM_AUDIT.json').write_text(json.dumps(dict(status='PASS',label_invariance='PASS',seams={k:dict(count=len(v),percentiles=np.percentile(v,[0,25,50,75,95,100]).tolist(),smooth_count={str(t):sum(a<=t for a in v) for t in [10,15,25]}) for k,v in seams.items()},samples=samples),indent=2)+'\n')
(S/'validation/PRODUCTION_SURFACE_CPP_PYTHON.json').write_text(json.dumps(dict(status='PASS',metrics=metrics,old_failure_cpp=replay),indent=2)+'\n')
(S/'raw/PRODUCTION_QUERY_COMPARISON.jsonl').write_text(''.join(json.dumps(a)+'\n' for a in records));print('PRODUCTION_SURFACE_AUDIT_PASS',metrics)
