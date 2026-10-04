from pathlib import Path
import json
import numpy as np


def write_vtk(folder,step,cloud,x,v,integrity,fields,face_traction):
    import pyvista as pv
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    poly=pv.PolyData(x)
    poly['reference_position_m']=cloud.X
    poly['displacement']=x-cloud.X;poly['velocity']=v
    poly['particle_volume_m3']=cloud.volume
    poly['FIXED_BASE']=cloud.fixed.astype(np.uint8)
    poly['EXPOSED_SURFACE']=cloud.exposed.astype(np.uint8)
    poly['INTERIOR']=(~cloud.exposed & ~cloud.fixed).astype(np.uint8)
    for key,value in fields.items():poly[key]=value
    poly.save(folder/f'particles_{step:04d}.vtp')
    bonds=pv.PolyData(x,lines=np.column_stack((np.full(len(cloud.pairs),2),cloud.pairs)).ravel())
    bonds.cell_data['integrity']=integrity
    bonds.cell_data['bond_damage']=1-integrity
    bonds.cell_data['broken']=(integrity==0).astype(np.uint8)
    bonds.save(folder/f'bonds_{step:04d}.vtp')
    surface=pv.PolyData(cloud.face_position)
    surface['solid_outward_normal']=cloud.face_normal
    surface['quadrature_area_m2']=cloud.face_area
    surface['particle_id']=cloud.face_particle
    surface['traction_Pa']=face_traction
    normal=np.sum(face_traction*cloud.face_normal,axis=1)[:,None]*cloud.face_normal
    surface['normal_traction_Pa']=normal
    surface['tangential_traction_Pa']=face_traction-normal
    surface.save(folder/f'traction_{step:04d}.vtp')


def write_collection(folder,steps):
    for kind in ['particles','bonds','traction']:
        text='<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>\n'
        text+='\n'.join(f'<DataSet timestep="{n}" group="" part="0" file="vtk/{kind}_{i:04d}.vtp"/>' for i,n in steps)
        (Path(folder)/f'{kind}.pvd').write_text(text+'\n</Collection></VTKFile>\n')


def write_json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
