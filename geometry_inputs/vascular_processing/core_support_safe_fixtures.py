"""Permanent deterministic seam, lever, bifurcation and wrap regressions."""
from pathlib import Path
import json,pickle,shutil
import numpy as np
import trimesh
from .support_inputs import save,ROOT
from .core_support_inputs import default_output
from .core_support_geometry import union,overlap,bed_roots
from .core_support_safe_geometry import thin_web_seam
from .core_support_slicing import evaluate
from .core_support_safe_removal import preserved_feature

NAMES=['continuous_roof_two_segments','staggered_roof_seam','thin_web_seam','long_lever_support','short_contact_piece','bifurcation_pull_bad','bifurcation_downward_remove_good','wrap_90','wrap_120','wrap_160','closed_ring_support']


def bifurcation_fixture():
    from .core_support_fixtures import fixture as original_fixture
    from .core_support_geometry import envelope,ruled
    vessel,branches=original_fixture('Y_bifurcation');offset=np.array([-115.,-110.,-6.])
    vessel.apply_translation(offset)
    for branch in branches:branch['xyz']=branch['xyz']+offset
    profile=dict(line_width=.22,support_line_width=.22,layer_height=.1,roof_gap_mm=.2)
    env=envelope(vessel,branches,profile,.2)
    support=union([ruled(b['xyz'][:,:2],b['roof_top'],.44,.44,.1)['mesh'] for b in env])
    return vessel,support


def sector(angle,radius=.9,length=2):
    # An annular prism, axis X. A 360-degree request has a genuine closed hole.
    if angle==360:
        m=trimesh.creation.annulus(r_min=radius,r_max=radius+.22,height=length,sections=96)
        m.apply_transform(trimesh.geometry.align_vectors([0,0,1],[1,0,0]));return m
    from shapely.geometry import Polygon
    theta=np.radians(np.linspace(-angle/2,angle/2,81));outer=np.c_[np.sin(theta)*(radius+.22),-np.cos(theta)*(radius+.22)];inner=np.c_[np.sin(theta)*radius,-np.cos(theta)*radius]
    m=trimesh.creation.extrude_polygon(Polygon(np.vstack([outer,inner[::-1]])),length,engine='triangle');m.apply_translation([0,0,-length/2]);m.apply_transform(trimesh.geometry.align_vectors([0,0,1],[1,0,0]));return m


def fixture(name):
    if name.startswith('wrap_'):return sector(float(name.split('_')[1]))
    if name=='closed_ring_support':return sector(360)
    if name.startswith('bifurcation_'):return bifurcation_fixture()[1]
    length=25 if name=='long_lever_support' else 2 if name=='short_contact_piece' else 12
    roof=trimesh.creation.box([length,.44,.3]);roof.apply_translation([0,0,5.15])
    if name=='short_contact_piece':return roof
    left=trimesh.creation.box([length/2-.22,.44,5]);left.apply_translation([-length/4-.11,0,2.5]);right=left.copy();right.apply_translation([length/2+.22,0,0]);m=union([roof,left,right])
    if name=='staggered_roof_seam':
        # Alternating shallow side notches, never a through-roof discontinuity.
        a=trimesh.creation.box([.22,.22,.1]);a.apply_translation([-.22,-.11,5.05]);b=a.copy();b.apply_translation([.44,.22,.1]);m=trimesh.boolean.difference([m,union([a,b])],engine='manifold')
    return m


def generate(out=None):
    out=Path(out or ROOT/'tests/data/vascular_core_support/safe_removal');out.mkdir(parents=True,exist_ok=True);rows=[]
    for name in NAMES:
        m=fixture(name);m.export(out/(name+'.stl'));rows.append(dict(name=name,watertight=bool(m.is_volume),path=str(out/(name+'.stl'))))
        if name.startswith('bifurcation_'):
            vessel,_=bifurcation_fixture();vessel.export(out/(name+'_vascular.stl'))
            rows[-1].update(vascular_path=str(out/(name+'_vascular.stl')),proposed_direction=[1,0,0] if name.endswith('bad') else [0,0,-1],expected_policy='REJECT_AXIAL_BIFURCATION_DRAG' if name.endswith('bad') else 'PREFER_DOWNWARD_RIGID_RELEASE')
    save(out/'manifest.json',dict(fixtures=rows,physical_validation=False));return rows


def actual_thin_web(out=None,width_lines=2):
    base=default_output();out=Path(out or base/'closure_v2/validation/actual_thin_web');out.mkdir(parents=True,exist_ok=True)
    source=base/'validation/actual_horizontal_tube';g=pickle.load((source/'geometry.pkl').open('rb'));p=g['profile'];r=json.loads((source/'geometry.json').read_text())
    # Add a real local curtain notch at x=115. Continuous roof unchanged.
    bottom=float(g['branches'][0]['roof_top'][0])-3*p['layer_height'];support,link=thin_web_seam(g['support'],[115,110,bottom],[1,0,0],p,bottom-.3,bottom,width_lines=width_lines)
    if width_lines>2:
        web=trimesh.creation.box([6*p['support_line_width'],width_lines*p['support_line_width'],3*p['layer_height']])
        web.apply_translation([115,110,bottom-1.5*p['layer_height']]);support=union([support,web])
    link.pop('cut_volume');g['support']=support
    for name in ['vascular_core.stl','vascular_print_transform.json','inverse_transform.json','notches.json']:shutil.copy2(source/name,out/name)
    support.export(out/'custom_support.stl');trimesh.util.concatenate([g['vascular'],support]).export(out/'vascular_core_plus_support_review.stl')
    with (out/'geometry.pkl').open('wb') as f:pickle.dump(g,f)
    r.update(candidate_id=out.name,folder=str(out),support_volume_mm3=float(support.volume),intersection_mm3=overlap(support,g['vascular']),bed_roots=bed_roots(support));save(out/'geometry.json',r)
    context=pickle.load((base/'source/context.pkl').open('rb'));e=evaluate(r,context['actual'],context['engine']);fp=pickle.load((out/'actual_footprints.pkl').open('rb'))['support']
    # Expected witness layers come from geometry/profile, never from the
    # surviving extrusion keys: a fully dropped intermediate layer must fail.
    levels=np.arange(p['first_layer_height'],support.bounds[1,2]+p['layer_height'],p['layer_height'])
    levels=levels[(levels>=bottom-.3+p['layer_height']/2)&(levels<=bottom+p['layer_height']/2)]
    link['witness_points']=[[115,110,float(z)] for z in levels]
    check=preserved_feature(link,fp,p)
    result=dict(weak_link=link,actual_preservation=check,evaluation=e,scope='SYNTHETIC_ONLY',vascular_source_unchanged=True,physical_validation=False)
    save(out/'thin_web_result.json',result);print('THIN_WEB_ACTUAL',check['passed'],e['coverage'],flush=True);return result

if __name__=='__main__':
    generate();actual_thin_web()
