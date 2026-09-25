"""Integrate only defined velocities; expose missing support rather than invent flow values."""
from pathlib import Path
import json,csv,sys
import numpy as np,yaml
from flow_geometry import FrozenSampler,Geometry,refine,quadrature
S=Path(__file__).resolve().parents[1]
cfg=yaml.safe_load((S/'configs/inlet_flux_audit.yaml').read_text())
flow=FrozenSampler(S/cfg['paths']['field']);geom=Geometry(S/cfg['paths']['geometry_arrays'])
m=json.loads((S/cfg['paths']['boundary_manifest']).read_text());inlet=next(p for p in m['ports'] if p['role']=='INLET_PORT')
normal=np.array(inlet['outward_unit_normal']);origin=np.array(inlet['center_m'])
offset=max(cfg['injection']['minimum_offset_voxels']*flow.dx,cfg['injection']['source_radius_support_upper_m']+flow.dx)
section=geom.section(origin-offset*normal,normal);np.savez_compressed(S/'geometry/INJECTION_SECTION.npz',**section)
tri=section['triangles'];records=[]
for level in range(max(cfg['quadrature']['refinement_levels'])+1):
    if level in cfg['quadrature']['refinement_levels']:
        points,weights=quadrature(tri);st,u=flow.query(points);valid=st==0;un=u@(-normal)
        qp=float(weights[valid]@np.maximum(un[valid],0));qb=float(weights[valid]@np.maximum(-un[valid],0));qn=float(weights[valid]@un[valid])
        areas={name:float(weights[st==code].sum()) for code,name in enumerate(['VALID','OUTSIDE','SOLID','MISSING','NONFINITE_POSITION'])}
        result={'level':level,'triangles':len(tri),'quadrature_points':len(points),'total_area_m2':float(weights.sum()),'status_area_m2':areas,'valid_area_fraction':areas['VALID']/weights.sum(),'partial_Q_positive_m3_s':qp,'partial_Q_net_m3_s':qn,'partial_Q_backflow_m3_s':qb,'partial_Q_backflow_fraction':qb/(qp+qb) if qp+qb else 0,'balance_residual_m3_s':qn-(qp-qb)}
        records.append(result);print(json.dumps(result),flush=True)
        if level==max(cfg['quadrature']['refinement_levels']):
            np.savez_compressed(S/'raw/INLET_QUADRATURE_FINE.npz',points=points,weights=weights,status=st,velocity=u,positive_normal_velocity=np.where(valid,np.maximum(un,0),np.nan))
            select=np.unique(np.r_[np.linspace(0,len(points)-1,min(2048,len(points)),dtype=int),np.flatnonzero(st!=0)[:100]])
            np.savetxt(S/'inputs/SAMPLER_COMPARISON_POINTS.txt',points[select],fmt='%.17g')
            np.savez_compressed(S/'raw/SAMPLER_PYTHON_REFERENCE.npz',points=points[select],status=st[select],velocity=u[select])
    if level<max(cfg['quadrature']['refinement_levels']):tri=refine(tri)
last,prev=records[-1],records[-2]
relative=abs(last['partial_Q_positive_m3_s']-prev['partial_Q_positive_m3_s'])/max(abs(last['partial_Q_positive_m3_s']),np.finfo(float).tiny)
coverage=last['valid_area_fraction'];complete=bool(coverage>=1-cfg['quadrature']['area_roundoff_tolerance_fraction'])
relative=float(relative);converged=bool(relative<=cfg['quadrature']['relative_Q_tolerance'])
reason=None if complete and converged else ('STOP_INLET_FLUX_INCOMPLETE_SAMPLER_COVERAGE' if not complete else 'STOP_INLET_FLUX_NOT_CONVERGED')
audit={'INLET_FLUX_INTEGRATION_STATUS':'PASS' if reason is None else 'BLOCKED','STOP_REASON':reason,'FLOW_PHYSICS_STATUS':'ENGINEERING_TRANSIENT_FIELD_ONLY','inward_unit_normal':(-normal).tolist(),'offset_m':offset,'offset_rule':'max(minimum_offset_voxels*dx, source_radius_support_upper_m+one_voxel_clearance)','clearance_m':flow.dx,'section_area_m2':float(section['area']),'quadrature':records,'highest_two_partial_Q_relative_difference':relative,'partial_integral_convergence_gate':converged,'full_section_velocity_coverage_gate':complete,'Q_positive_m3_s':last['partial_Q_positive_m3_s'] if complete else None,'Q_net_m3_s':last['partial_Q_net_m3_s'] if complete else None,'Q_backflow_m3_s':last['partial_Q_backflow_m3_s'] if complete else None,'Q_backflow_fraction':last['partial_Q_backflow_fraction'] if complete else None,'partial_integrals_are_not_full_inlet_flux':not complete,'missing_velocity_policy':'No renormalization, no extrapolation, no physical zero assigned to invalid query. Status SOLID means at least one interpolation corner is not fluid; it does not prove quadrature point is outside the continuous lumen.'}
(S/'validation/INLET_FLUX_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps({k:v for k,v in audit.items() if k!='quadrature'},indent=2))
