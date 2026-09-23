#!/usr/bin/env python3
"""Independent extended-precision check of the rejected area invariant."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import pyvista as pv
from fem3d.audit import write_json,timestamp
b=json.loads((ROOT/'inputs/stage01_5/baseline_contract.json').read_text())
source=np.load(ROOT/'inputs/stage01/tagged_surface_si.npz')
def measure(data,port):
    p=np.asarray(data['points_m'],dtype=np.longdouble);t=data['triangles'][data['facet_tags']==port['surface_entity_id']]
    c=p[t];v=np.cross(c[:,1]-c[:,0],c[:,2]-c[:,0])/2
    area=np.sqrt(np.sum(v*v,axis=1)).sum(dtype=np.longdouble)
    n=np.asarray(port['outward_normal'],dtype=np.longdouble);n/=np.sqrt(n@n)
    vector=v.sum(axis=0,dtype=np.longdouble)
    projected=vector@n
    return area,projected,vector
records={}
for port in b['ports']:
    name=port['name'];old,oldproj,oldvec=measure(source,port)
    row={'source_area_m2':float(old),'source_projected_area_m2':float(oldproj),'source_scalar_minus_projected_relative':float((old-oldproj)/old),'source_max_plane_deviation_m':port['max_plane_deviation_m'],'candidates':{}}
    for candidate in ('candidate_A','candidate_B','candidate_C'):
        data=np.load(ROOT/'outputs/stage01_5'/candidate/'surface/tagged_surface_si.npz')
        new,projected,vec=measure(data,port)
        row['candidates'][candidate]={'extended_precision_relative_error_to_contract':float(abs(new-np.longdouble(port['area_m2']))/np.longdouble(port['area_m2'])),'signed_relative_area_change_vs_source':float((new-old)/old),'projected_area_relative_change':float((projected-oldproj)/oldproj),'vector_area_relative_change':float(np.sqrt(np.sum((vec-oldvec)**2))/oldproj),'scalar_minus_projected_relative':float((new-projected)/projected)}
        assert row['candidates'][candidate]['extended_precision_relative_error_to_contract']>1e-12
    records[name]=row
write_json(ROOT/'reports/stage01_5/area_invariant_analysis.json',{'timestamp':timestamp(),'status':'CONFIRMED GEOMETRY GATE FAILURE','arithmetic':{'dtype':'numpy.longdouble','mantissa_bits':int(np.finfo(np.longdouble).nmant),'epsilon':float(np.finfo(np.longdouble).eps)},'source_vtp_coordinate_dtype':str(pv.read(ROOT/'inputs/stage00/source_surface_um.vtp').points.dtype),'interpretation':'Fixed noncoplanar rim preserves oriented vector area and projected polygon area. Scalar 3D triangle area depends on triangulation. Interior points on the fixed plane change this scalar area by more than the frozen tolerance. This is not a relaxable floating-point sum discrepancy.','ports':records})
print('All 12 port/candidate area failures independently confirmed in extended precision')
