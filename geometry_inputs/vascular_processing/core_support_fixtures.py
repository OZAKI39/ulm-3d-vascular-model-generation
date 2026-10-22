"""Deterministic vascular-only synthetic geometry; no casting parts."""
from pathlib import Path
import numpy as np
from .support_rfus_fixtures import tube,branch
from .core_support_geometry import union,generate
from .core_support_slicing import evaluate
from .support_inputs import save

NAMES=['horizontal_tube','curved_tube','S_tube','Y_bifurcation','crossing_tubes','stacked_tubes','four_port_tree','near_vertical_tube','long_horizontal_port']


def fixture(name):
    if name not in NAMES:raise ValueError(name)
    x=np.linspace(100,130,151);y=np.full(len(x),110.);z=np.full(len(x),10.)
    if name=='curved_tube':y+=3*np.sin(np.linspace(0,np.pi,len(x)))
    if name=='S_tube':y+=2*np.sin(np.linspace(0,2*np.pi,len(x)))
    if name=='near_vertical_tube':x=np.linspace(110,111,len(x));z=np.linspace(5,25,len(x))
    if name=='long_horizontal_port':x=np.linspace(90,140,len(x))
    xyz=np.c_[x,y,z];bs=[branch(xyz,name=='long_horizontal_port','MAIN')]
    if name in ['Y_bifurcation','four_port_tree']:
        daughter=np.c_[np.linspace(115,130,81),np.linspace(110,122,81),np.full(81,10.)]
        bs=[branch(xyz[:76],False,'PARENT'),branch(xyz[75:],name=='four_port_tree','DAUGHTER_1'),branch(daughter,name=='four_port_tree','DAUGHTER_2')]
        if name=='four_port_tree':
            bs.extend([branch([[100,110,10],[92,110,10]],True,'PORT_I1'),branch([[115,110,10],[115,100,10]],True,'PORT_O3')]);bs[0]['is_port']=False
    if name=='crossing_tubes':bs.append(branch([[115,102,5],[115,118,5]],False,'LOWER'))
    if name=='stacked_tubes':bs.append(branch([[100,110,5],[130,110,5]],False,'LOWER'))
    mesh=union([tube(b['xyz']) for b in bs]);return mesh,bs


def generate_all(out,profile):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    for name in NAMES:
        mesh,bs=fixture(name);mesh.export(out/(name+'.stl'));save(out/(name+'.json'),dict(name=name,branches=bs,vascular_triangle_count=len(mesh.faces),scale_factor=1.));rows.append(dict(name=name,watertight=bool(mesh.is_volume),file=str(out/(name+'.stl'))))
    save(out/'manifest.json',dict(fixtures=rows,profile=profile,physical_validation=False));return rows


def actual_horizontal(out,profile,actual,engine):
    mesh,bs=fixture('horizontal_tube');row=dict(orientation_id='SYNTHETIC_IDENTITY',transform=np.eye(4),rotation=np.eye(3))
    r=generate(mesh,bs,row,profile,Path(out)/'actual_horizontal_tube',gap=.2,lines=2,roof_lines=2)
    # A perfectly horizontal isolated tube leaves globally empty layers across
    # its intentional gap. Bambu rejects that fixture without an independent
    # layer witness. The witness is outside all tube/support target projections
    # and receives no underbody coverage credit.
    import pickle,trimesh
    folder=Path(r['folder']);g=pickle.load((folder/'geometry.pkl').open('rb'))
    witness=trimesh.creation.box([1.,1.,12.]);witness.apply_translation([136,118,6.])
    g['support']=union([g['support'],witness]);g['support'].export(folder/'custom_support.stl');trimesh.util.concatenate([g['vascular'],g['support']]).export(folder/'vascular_core_plus_support_review.stl')
    with (folder/'geometry.pkl').open('wb') as f:pickle.dump(g,f)
    e=evaluate(r,actual,engine);save(Path(out)/'horizontal_actual_result.json',dict(length_mm=30.,layer_witness=dict(center=[136,118,6],size_mm=[1,1,12],credited_as_underbody=False,scope='SYNTHETIC_FIXTURE_ONLY'),**e));return r,e
