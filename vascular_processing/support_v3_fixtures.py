"""Permanent synthetic geometry for s1-8 V3 regressions (no physical claims)."""
from pathlib import Path
import numpy as np
import trimesh
from .support_v3_geometry import curtain, union
from .support_v2_primitives import loft
from .support_inputs import save

NAMES=['straight_tube_long_rib','curved_tube_segmented_rib','S_tube_rib','bifurcation_rib',
       'external_port_rib','rib_window','rib_roof_transition','rib_cut_notch','rib_tool_access',
       'rib_fragment_removal','rib_self_support','one_line_rib','two_line_rib']


def fixture(name,profile):
    if name not in NAMES:raise ValueError(name)
    x=np.linspace(0,30,121);y=np.zeros_like(x)
    if name=='curved_tube_segmented_rib':y=3*np.sin(x/30*np.pi)
    if name=='S_tube_rib':y=2*np.sin(x/30*2*np.pi)
    lines=1 if name=='one_line_rib' else 2
    rib=curtain(np.c_[x,y],np.full(len(x),10.),0,profile,lines,6*profile['line_width'],5.,3.)
    obstacles=[]
    if name=='bifurcation_rib':
        side=curtain(np.c_[np.linspace(12,20,41),np.linspace(1,9,41)],np.full(41,10.),0,profile,2,6*profile['line_width'],3.,3.)
        rib['mesh']=union([rib['mesh'],side['mesh']])
    if name=='external_port_rib':
        obstacles=[loft([-4,0,0],[-4,0,12],[2,12])]
    return rib,obstacles


def generate(out,profile):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    for name in NAMES:
        rib,obstacles=fixture(name,profile);rib['mesh'].export(out/f'{name}.stl')
        if obstacles:union(obstacles).export(out/f'{name}_obstacles.stl')
        rows.append(dict(fixture=name,volume_mm3=float(rib['mesh'].volume),windows=len(rib['windows']),
                         cut_points=len(rib['cuts']),line_width_mm=profile['line_width'],physical_validation=False))
    save(out/'manifest.json',dict(architecture='FDM_SEGMENTED_RIB_SUPPORT_V3',fixtures=rows))
    return rows
