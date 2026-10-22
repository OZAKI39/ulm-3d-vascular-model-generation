"""Refine shared trunks into designed, printable removal sections."""
import copy,pickle,shutil
from pathlib import Path
import numpy as np
import trimesh
from . import support_v2_primitives as p
from .support_v2_pipeline import plan_candidate
from .support_inputs import save

def segmented_pillar(a,b,width,profile):
    a=np.array(a,float);b=np.array(b,float);h=profile['layer_height'];w=profile['line_width']
    pitch=p.quantize(3*profile['bridge_limit_mm'],h);neck_width=p.quantize(max(3*w,width*.5),w)
    flare_h=p.quantize((width-neck_width)/2/.6,h);neck_h=3*h;shapes=[];cuts=[];records=[]
    def point(z):return a+(b-a)*(z-a[2])/(b[2]-a[2])
    def add(mesh,kind,dimension):shapes.append(mesh);records.append(dict(kind=kind,width_mm=dimension))
    previous=a
    for z in np.arange(np.ceil((a[2]+neck_h+2*h)/pitch)*pitch,b[2]-flare_h-neck_h,pitch):
        bottom=point(z-neck_h/2);top=point(z+neck_h/2);flaretop=point(z+neck_h/2+flare_h)
        if bottom[2]-previous[2]<h:continue
        add(p.pillar(previous,bottom,width),'VERTICAL_PILLAR',width)
        add(p.loft(bottom,top,neck_width),'BREAKAWAY_SEGMENT_NECK',neck_width)
        add(p.buttress(top,flaretop,neck_width,width),'BUTTRESS',width)
        cuts.append(dict(xyz=point(z).tolist(),width_mm=neck_width+.08,height_mm=.06,role='TRUNK_BREAK',neck_width_mm=neck_width,neck_height_mm=neck_h))
        previous=flaretop
    add(p.pillar(previous,b,width),'VERTICAL_PILLAR',width)
    return shapes,records,cuts

def refine(record,inputs,profile,cfg,mold,v,box,branches,topology,out):
    folder=Path(record['folder']);cache=folder/'scaffold_local_cache.pkl'
    with cache.open('rb') as f:graph,roofs,pieces,approved=pickle.load(f)
    if graph.get('trunk_segmentation_refinement_completed'):return record
    archive=folder/'before_trunk_segmentation';archive.mkdir(exist_ok=True)
    for name in ['scaffold_local_cache.pkl','custom_support.stl','mold_with_custom_support.stl','candidate.json','support_graph_v2.json']:
        if (folder/name).exists():shutil.copy2(folder/name,archive/name)
    rebuilt=[];new_records=[];extra=[]
    for mesh,meta in zip(pieces,graph['primitives']):
        if meta['kind']=='VERTICAL_PILLAR' and 'root_id' in meta:
            lo,hi=mesh.bounds;bottom=mesh.vertices[np.abs(mesh.vertices[:,2]-lo[2])<1e-6].mean(0);top=mesh.vertices[np.abs(mesh.vertices[:,2]-hi[2])<1e-6].mean(0)
            if top[2]-bottom[2]>.5:
                shapes,records,cuts=segmented_pillar(bottom,top,meta['width_mm'],profile)
                for shape,new in zip(shapes,records):rebuilt.append(shape);new_records.append(dict(**new,root_id=meta['root_id'],supported_load=meta.get('supported_load',1)))
                extra.extend({**c,'root_id':meta['root_id']} for c in cuts);continue
        rebuilt.append(mesh);new_records.append({k:v for k,v in meta.items() if k not in ['primitive_id','bbox_mm','volume_mm3']})
    for i,(mesh,meta) in enumerate(zip(rebuilt,new_records)):meta.update(primitive_id=f'P{i:04}',bbox_mm=mesh.bounds.tolist(),volume_mm3=float(mesh.volume))
    graph['primitives']=new_records
    for cut in extra:cut['neck_id']=f"N{len(graph['segment_necks']):04}";graph['segment_necks'].append(cut)
    graph['trunk_segmentation_refinement_completed']=True;graph['added_trunk_break_necks']=len(extra)
    staging=folder/'refinement_work';work=staging/'candidates'/record['candidate_id'];work.mkdir(parents=True,exist_ok=True)
    with (work/'scaffold_local_cache.pkl').open('wb') as f:pickle.dump((graph,roofs,rebuilt,approved),f)
    # Generate and audit a complete replacement before touching the candidate.
    # A failed precision/Boolean export must not mix old STL with a new graph.
    r=plan_candidate(staging,record['orientation'],inputs,profile,cfg,mold,v,box,branches,topology,True)
    r['refinement']=dict(iteration=3,added_trunk_break_necks=len(extra),previous_support_volume_mm3=record['support_volume_mm3'],contact_positions_unchanged=True,orientation_unchanged=True,reason='Shared dense trunks may otherwise remain mechanically connected after root/head neck cuts; insert common-height designed sections')
    for path in work.iterdir():
        if path.is_dir():shutil.copytree(path,folder/path.name,dirs_exist_ok=True)
        else:shutil.copy2(path,folder/path.name)
    r['folder']=str(folder);save(folder/'candidate.json',r);return r
