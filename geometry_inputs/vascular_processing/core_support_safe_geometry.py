"""Local cuttable seams for s1-11; gated by frozen real coverage.

The roof itself stays continuous during printing. The smaller connection is
in the curtain below the roof. Its cross-section is only a geometric proxy.
"""
from __future__ import annotations
import copy,pickle,json,shutil
from pathlib import Path
import numpy as np
import trimesh
from .core_support_geometry import boolean,union,overlap,bed_roots
from .support_v2_mesh import support_stl
from .core_support_slicing import evaluate
from .core_support_closure import closure_pass
from .core_support_safe_removal import preserved_feature,cross_section_area,weak_ratio
from .support_inputs import save,table,digest


def thin_web_seam(support,center,tangent,profile,body_top,roof_bottom,width_lines=2):
    """Notch the curtain body only; retain a roof-side connecting web."""
    w=profile['support_line_width'];h=profile['layer_height'];t=np.asarray(tangent,dtype=float);t/=np.linalg.norm(t)
    M=np.eye(4);M[:2,0]=t[:2];M[:2,1]=[-t[1],t[0]];M[:3,3]=[*np.asarray(center)[:2],(body_top-h)/2]
    cutter=trimesh.creation.box([2*w,8*w,body_top+h],transform=M)
    result=boolean([support,cutter],'difference')
    return result,dict(center=np.asarray(center).tolist(),weak_link_width_mm=width_lines*w,weak_link_thickness_mm=roof_bottom-body_top,weak_link_area_mm2=width_lines*w*(roof_bottom-body_top),quantization=dict(width_lines=width_lines,thickness_layers=(roof_bottom-body_top)/h),cut_volume=cutter)


def removal_preserves_coverage(before,after):
    return closure_pass(before) and closure_pass(after)


def optimize_after_closure(record,before,out,actual,engine):
    """One bounded seam trial, followed by real slicing and strict rollback.

    Called only after real coverage closure. Never changes source geometry and
    never promotes a removal trial that loses coverage or a critical feature.
    """
    if not closure_pass(before):raise ValueError('REMOVAL_GEOMETRY_LOCKED_UNTIL_100_PERCENT_COVERAGE')
    source=Path(record['folder']);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    g=pickle.load((source/'geometry.pkl').open('rb'));support=g['support'];v=g['vascular'];p=g['profile'];h=p['layer_height'];w=p['support_line_width'];links=[]
    samples=np.vstack([b['xyz'] for b in g['branches']]);tangents=np.vstack([b['tangent'] for b in g['branches']])
    from scipy.spatial import cKDTree
    for n in g['notches']:
        roof_bottom=n['plane_origin'][2]+.15;top=roof_bottom-3*h
        trial,link=thin_web_seam(support,n['point'],n['tangent'],p,top,roof_bottom)
        cutter=link.pop('cut_volume')
        if overlap(cutter,v)>=1e-7:continue
        support=trial;idx=cKDTree(samples).query(n['point'])[1]
        area=cross_section_area(v,samples[idx],tangents[idx]);link.update(link_id=n['notch_id'],nearest_vascular_area_mm2=area,geometric_weak_link_ratio=weak_ratio(link['weak_link_width_mm'],link['weak_link_thickness_mm'],area))
        first=p['first_layer_height'];levels=np.arange(first,support.bounds[1,2]+h,h);levels=levels[(levels>=top+h/2)&(levels<=roof_bottom+h/2)]
        link['witness_points']=[[*np.asarray(n['point'])[:2],float(z)] for z in levels];links.append(link)
    support=support_stl(support);g['support']=support;g['weak_links']=links
    for name in ['vascular_core.stl','vascular_print_transform.json','inverse_transform.json','notches.json']:shutil.copy2(source/name,out/name)
    support.export(out/'custom_support.stl');trimesh.util.concatenate([v,support]).export(out/'vascular_core_plus_support_review.stl')
    shifted=support.copy();shifted.apply_translation([0,0,-3*h]);g['roofs']=boolean([support,shifted],'difference');g['curtains']=boolean([support,g['roofs']],'difference')
    for name in ['roofs','curtains']:g[name].export(out/('support_'+name+'.stl'))
    with (out/'geometry.pkl').open('wb') as f:pickle.dump(g,f)
    r=copy.deepcopy(record);r.update(candidate_id=out.name,folder=str(out),support_volume_mm3=float(support.volume),intersection_mm3=overlap(support,v),bed_roots=bed_roots(support));save(out/'geometry.json',r)
    e=evaluate(r,actual,engine);fp=pickle.load((out/'actual_footprints.pkl').open('rb'))['support']
    for link in links:link['actual_preservation']=preserved_feature(link,fp,p)
    g['weak_links']=links
    with (out/'geometry.pkl').open('wb') as f:pickle.dump(g,f)
    accepted=removal_preserves_coverage(before,e) and bool(links) and all(x['actual_preservation']['passed'] for x in links) and r['intersection_mm3']<1e-7
    save(out/'weak_link_preservation.json',dict(links=links,accepted=accepted,coverage_preserved=removal_preserves_coverage(before,e),source_sha256=digest(out/'vascular_core.stl'),physical_validated=False))
    return (r,e) if accepted else (record,before)
