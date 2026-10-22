"""Deterministic tolerance alternatives; engineering stress cases, not calibration."""
import copy,itertools
from pathlib import Path
import numpy as np
from .support_inputs import save,table
from .casting_fixture_audit import perturbation_fit
from .casting_fixture_geometry import source_layout,build_frame,export_shape,mesh_qc,shape_mesh


def evaluate(inputs,cfg,line_width,out,export_alternatives=True):
    out=Path(out);layout=source_layout(inputs,cfg,line_width);ports=layout['ports'];secondary=[i for i,p in enumerate(ports) if p['role']!='MASTER_DATUM'];errors=[]
    radial=.15;axial=.75
    for signs in itertools.product((-1,0,1),repeat=len(secondary)):
        e=np.zeros((len(ports),3))
        for i,sign in zip(secondary,signs):e[i]=np.array(ports[i]['axis'])*sign*axial
        errors.append(('AXIAL_'+str(signs),e))
    for signs in itertools.product((-1,1),repeat=len(ports)):
        e=np.array([[0.,0,sign*radial] for sign in signs]);errors.append(('RADIAL_'+str(signs),e))
        combined=e.copy()
        for i in secondary:combined[i]+=np.array(ports[i]['axis'])*signs[i]*axial
        errors.append(('COMBINED_'+str(signs),combined))
    rows=[];details=[]
    for multiplier in cfg['port_radial_clearance_line_multipliers']:
        for floating in cfg['axial_float_candidates_mm']:
            variant=copy.deepcopy(ports)
            for p in variant:p['radial_clearance_mm']=line_width*multiplier;p['axial_float_mm']=0 if p['role']=='MASTER_DATUM' else floating
            success=0;fixed=0
            for name,e in errors:
                good=perturbation_fit(variant,e,True);bad=perturbation_fit(variant,e,False);success+=good;fixed+=bad;details.append(dict(radial_clearance_mm=line_width*multiplier,axial_float_mm=floating,case=name,floating_pass=good,all_fixed_pass=bad,errors_mm=e.tolist()))
            rows.append(dict(radial_clearance_mm=line_width*multiplier,axial_float_mm=floating,test_count=len(errors),floating_pass_count=success,fully_fixed_pass_count=fixed,all_stress_cases_passed=success==len(errors)))
    winners=[r for r in rows if r['all_stress_cases_passed']];selected=min(winners,key=lambda r:(r['radial_clearance_mm'],r['axial_float_mm'])) if winners else None
    result=dict(physical_tolerance_calibrated=False,radial_stress_mm=radial,axial_stress_mm=axial,assumptions='Explicit deterministic engineering stress envelope; not measured PLA/ABS error',no_monte_carlo=True,candidate_count=len(rows),rows=rows,smallest_passing_candidate=selected,current_candidate=dict(radial_clearance_mm=line_width*cfg['selected_radial_clearance_line_multiplier'],axial_float_mm=cfg['selected_axial_float_mm']))
    table(out/'tolerance_candidates.csv',rows);table(out/'deterministic_perturbations.csv',details);save(out/'tolerance_summary.json',result)
    if export_alternatives:
        for multiplier in cfg['port_radial_clearance_line_multipliers']:
            if multiplier==cfg['selected_radial_clearance_line_multiplier']:continue
            alternative=copy.deepcopy(cfg);alternative['selected_radial_clearance_line_multiplier']=multiplier;l=source_layout(inputs,alternative,line_width);g=build_frame(l,alternative,line_width);folder=out/'geometry_candidates'/('radial_'+format(line_width*multiplier,'.2f'))
            export_shape(g['frame'],folder/'frame',cfg)
            for k,v in g['closures'].items():export_shape(v,folder/f'closure_{k}',cfg)
            save(folder/'parameters.json',dict(radial_clearance_mm=line_width*multiplier,axial_float_mm=cfg['selected_axial_float_mm'],role='CLEARANCE_REVIEW_CANDIDATE_NOT_SELECTED_PRINT_FILES'))
    return result
