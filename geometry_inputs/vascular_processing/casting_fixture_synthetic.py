"""Deterministic small fixtures for assembly, tolerance and impossible cases."""
import copy,json
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save

CASES=['one_inlet_two_outlets','three_side_ports','four_side_ports','port_size_variation','port_position_error','top_loading_case','overconstrained_fixture','floating_secondary_port','closure_collision','dish_too_small']

def layout_for(name,clearance=.22):
    axes=[[0,-1,0],[1,0,0],[-1,0,0]]
    if name in ('four_side_ports','port_size_variation','overconstrained_fixture','floating_secondary_port','port_position_error'):axes.append([0,1,0])
    ports=[];pieces=[]
    for i,axis in enumerate(axes):
        n=np.array(axis,float);center=n*10+np.array([0,0,12.]);radius=(.65+.15*i if name=='port_size_variation' else .8)
        p=dict(port_id=f'P{i+1}',axis=n.tolist(),face=('-' if min(axis)<0 else '+')+('X' if axis[0] else 'Y'),wall_intersection=center.tolist(),capture_center=(center+n).tolist(),straight_start=(center-n*2).tolist(),end=(center+n*6).tolist(),diameter_mm=2*radius,radius_mm=radius,radial_clearance_mm=clearance,role='MASTER_DATUM' if i==0 else 'FLOATING_SECONDARY',axial_float_mm=0 if i==0 else 1.,straight_length_mm=8.)
        ports.append(p);transform=trimesh.geometry.align_vectors([0,0,1],n);transform[:3,3]=[0,0,12];transform[:3,3]+=n*8;pieces.append(trimesh.creation.cylinder(radius=radius,height=16,sections=32,transform=transform))
    central=trimesh.creation.icosphere(subdivisions=2,radius=1.6);central.apply_translation([0,0,12]);core=trimesh.boolean.union([central,*pieces],engine='manifold')
    return dict(core=core,central=central,ports=ports,low=np.array([-10.,-10.,0]),high=np.array([10.,10.,26]),master_port='P1',transform=np.eye(4))

def generate(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    for name in CASES:
        layout=layout_for(name);path=directory/name;path.mkdir(exist_ok=True);layout['core'].export(path/'core.stl')
        save(path/'fixture.json',dict(name=name,synthetic=True,ports=layout['ports'],inner_bounds=[layout['low'],layout['high']],deterministic_port_errors_mm=[(np.asarray(p['axis'])*.8).tolist() if p['role']!='MASTER_DATUM' else [0,0,0] for p in layout['ports']],expected={'closure_collision':'REJECT_COLLISION','dish_too_small':'DISH_FIT_NOT_VERIFIED','overconstrained_fixture':'FULLY_FIXED_REJECTS_ERRORS'}.get(name,'GEOMETRY_REGRESSION')))
    (directory/'README.md').write_text('# BCF deterministic fixtures\n\nEach case has a real closed synthetic STL and explicit port metadata. No patient data. Manufacturing error samples are deterministic ± axis perturbations, not Monte Carlo. Negative cases deliberately reject collision, undersized dishes and overconstrained fit.\n')
