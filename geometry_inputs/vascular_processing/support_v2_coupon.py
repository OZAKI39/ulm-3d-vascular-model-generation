"""Gap, bridge, rail, Y and fracture experiments; no physical validation claim."""
import json,zipfile
import numpy as np
import trimesh
from . import support_v2_primitives as p
from .support_tree_planner import strut
from .support_inputs import save,table
from .support_custom_qc import canonical_stl,feature_preservation
from .support_slicer_adapter import discover_engine,engine_profiles
from .vascular_support_pipeline import sliced
from .support_qc import parse_actual_paths
from .bambu_casting_mold import sample_segments

def create(out,profile):
    out.mkdir(parents=True,exist_ok=True);w=profile['line_width'];h=profile['layer_height'];items=[];rows=[];features=[]
    specs=[('ROOF_GAP_NECK',gap,n,None) for gap in [.10,.15,.20,.25] for n in [2,3,4]]+ [('BRIDGE',.2,3,s) for s in [2,4,6,8]]+[('RAIL',.2,3,2),('Y_MERGE',.2,3,None)]
    for i,(kind,gap,lines,span) in enumerate(specs):
        xy=np.array([35+(i%6)*22,35+(i//6)*22]);base=p.loft([0,0,0],[0,0,2],[18,18]);parts=[base];targets=[]
        if kind=='ROOF_GAP_NECK':
            top=9.4-gap;roof=p.roof_pad([0,0,9.4],8*w,4*w,3*h,gap)[0];neckw=lines*w;flareh=p.quantize((8*w-neckw)/2/.6,h);flaretop=np.array([0,0,top-3*h]);flarebottom=flaretop-[0,0,flareh];neckcenter=flarebottom-[0,0,1.5*h]
            stem=p.pillar([0,0,2],[0,0,neckcenter[2]-1.5*h],8*w)
            parts.extend([stem,p.neck(neckcenter,neckw,3*h),p.buttress(flarebottom,flaretop,neckw,[8*w,4*w]),roof])
            parts.extend([p.pillar([x,0,2],[x,0,10],8*w) for x in [-4,4]])
            parts.append(strut([-4.5,0,10],[4.5,0,10],1.2))
            targets=[('ROOF_PAD',[0,0,top-h/2],4*w),('BREAKAWAY_SEGMENT_NECK',neckcenter,neckw),('SMALL_ABS_TUBE',[0,0,10],1.2)]
        elif kind in ['BRIDGE','RAIL']:
            width=8*w;distance=span+width;a=np.array([-distance/2,0,7.]);b=np.array([distance/2,0,7.])
            parts.extend([p.pillar([q[0],0,2],q,width) for q in [a,b]]);parts.append(p.horizontal_member(a,b,width,4*h))
            targets=[('SHORT_BRIDGE' if kind=='BRIDGE' else 'LONGITUDINAL_RAIL',[0,0,7],width)]
            if kind=='RAIL':
                for x in [-1.,0.,1.]:parts.append(p.pillar([x,0,7],[x,0,8],4*w))
                parts.append(strut([-3,0,8.8],[3,0,8.8],1.2));parts.extend([p.pillar([x,0,2],[x,0,8.8],8*w) for x in [-3,3]])
        else:
            parts.append(p.y_merge([-2,0,10],[2,0,10],[0,0,5],[0,0,2],6*w));targets=[('Y_MERGE',[0,0,5],6*w),('SELF_SUPPORTING_BRANCH',[-1,0,7.5],6*w)]
            for x in [-2,2]:
                parts.append(p.roof_pad([x,0,10.3],6*w,6*w,3*h,0)[0]);targets.append(('ROOF_PAD',[x,0,10.25],6*w))
        mesh=canonical_stl(trimesh.boolean.union(parts,engine='manifold'));mesh.apply_translation([*xy,0]);items.append(mesh)
        for name,xyz,width in targets:features.append(dict(feature_id=f'C{i:02}_{name}_{len(features)}',kind=name,xyz=(np.asarray(xyz)+[*xy,0]).tolist(),width_mm=width))
        rows.append(dict(coupon_id=f'C{i:02}',kind=kind,center_xy_mm=xy.tolist(),roof_gap_mm=gap,neck_width_mm=lines*w,bridge_free_span_mm=span,physical_test_completed=False))
    mesh=canonical_stl(trimesh.util.concatenate(items));mesh.export(out/'fdm_support_coupon_v2.stl');table(out/'coupon_layout.csv',rows);save(out/'coupon_features.json',features)
    save(out/'coupon_manifest.json',dict(cases=rows,bridge_limit_status=profile['bridge_calibration_status'],gap_status=profile['gap_status'],physical_test_completed=False))
    (out/'coupon_instructions.md').write_text('# FDM V2 试片\n\n18 个单元包括 4 种顶垫间隙 × 3 种分段颈宽，以及 2/4/6/8 mm 桥、短轨道、Y 合并和小直径 ABS 圆杆。位置见 coupon_layout.csv。使用 X1C 0.2 mm、ABS Vessel、0.10 mm 层高，关闭自动支撑。\n\n记录桥的下垂、各细颈是否完整、实际断裂位置、移除力和圆杆损伤；每种参数至少重复三次。桥跨 4/6/8 mm 是敏感性试验，不代表已经通过材料校准。0.15/0.25 mm 几何间隙在 0.10 mm 层高下会量化，按真实刀路确认实际间隙。只有实物试验通过并填写匹配工艺的 config/fdm_bridge_calibration.json，程序才会采用物理校准的桥接上限。\n')
    return mesh,rows,features

def run_coupon(out,profile,discovery):
    mesh,rows,features=create(out,profile);engine=discover_engine('Bambu',discovery,out);profiles=engine_profiles(discovery,'Bambu',{'enable_support':'0'})
    source=out/'fdm_support_coupon_v2.stl';r=sliced(engine,profiles,source,source,out/'bambu',900)
    if r.get('slicer_success'):
        paths=parse_actual_paths(r['archive']);preserve=feature_preservation(sample_segments(paths['segments'][paths['kind']==0],.1),features,profile)
        with zipfile.ZipFile(r['archive']) as archive:settings=json.loads(archive.read('Metadata/project_settings.config'))
        r['feature_preservation']=preserve;r['auto_support_off']=str(settings.get('enable_support')) in ['0','false','False'];r['native_support_toolpaths']=int((paths['kind']==1).sum());table(out/'feature_preservation.csv',preserve['rows'])
        from .support_v2_validation import bridge_paths,bridge_preservation
        graph=[]
        for cell in rows:
            if cell['kind'] not in ['BRIDGE','RAIL']:continue
            x,y=cell['center_xy_mm'];width=8*profile['line_width'];half=(cell['bridge_free_span_mm']+2*width)/2
            graph.append(dict(primitive_id=cell['coupon_id'],kind='SHORT_BRIDGE',bbox_mm=[[x-half,y-width/2,6.8],[x+half,y+width/2,7.2]],bridge_direction_xy=[1.,0.],span_mm=cell['bridge_free_span_mm']))
        r['bridge_toolpaths']=bridge_preservation(bridge_paths(r['archive']),{'primitives':graph},paths['segments'][paths['kind']==0],profile);table(out/'bridge_toolpath_preservation.csv',r['bridge_toolpaths']['rows'])
    save(out/'bambu_validation.json',r);print('V2 COUPON',r.get('slicer_success'),flush=True);return r
