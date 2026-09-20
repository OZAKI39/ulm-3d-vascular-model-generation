#!/usr/bin/env python3
"""Gated PETSc comparison and frozen time policy; no parameter tuning."""
import argparse, difflib, json, shutil, subprocess, sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import REPORT,load,compare_production_xml,validate_petsc_xml,linear_gate,nonlinear_gate,parse_boundary,write_csv,require_frozen_inputs
from sv_validation.sv11_runtime import OPTIONS,petsc_ls,run_solver
from sv_validation.provenance import sha256,write_json
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.validation import integral_agreement,steady_state

def measure(case,name,run):
    policy=json.loads((ROOT/'configs/time_policy.json').read_text());baseline=load('flow_qc','sv1')
    m=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',baseline['Q_target_m3_s'],policy['Umean_m_s'])
    native=case/'4-procs/B_NS_Velocity_flux.txt'
    boundary=parse_boundary(native,m.Q) if native.exists() else []
    states=[];last=None;intervals=[]
    for path in sorted(case.glob('4-procs/result_*.vtu')):
        u,p=m.read(path);q=m.measure(u,p);q.pop('outlet_fractions',None)
        step=int(path.stem.split('_')[-1]);q.update(step=step,path=str(path.relative_to(ROOT)),sha256=sha256(path),source='actual_vtu_surface_integration')
        if last:
            intervals.append({'step':step,'E_u':m.velocity_l2(u-last[0])/(m.velocity_l2(u)+np.finfo(float).tiny),
                              'E_Q':max(abs(q['signed_outward_boundary_flows_m3_s'][r]-last[1]['signed_outward_boundary_flows_m3_s'][r])/m.Q for r in q['signed_outward_boundary_flows_m3_s'])})
        match=next((b for b in boundary if b['step']==step),None)
        if match:
            flows={r:match['signed_outward_flows_m3_s'][r] for r in q['signed_outward_boundary_flows_m3_s']}
            q['native_flux_difference_over_Q']=integral_agreement(flows,q['signed_outward_boundary_flows_m3_s'],m.Q)
        states.append(q);last=(u,q)
    result={'states':states,'boundary_history':boundary,'steady_intervals':intervals,
            'steady_reached':steady_state([r['E_u'] for r in intervals],[r['E_Q'] for r in intervals]) if len(intervals)>=5 else False,
            'completed_steps':max((r['step'] for r in run['history']['linear_solves']),default=0)}
    write_json(REPORT/(name+'_qc.json'),result)
    if boundary:write_csv(REPORT/(name+'_boundary.csv'),boundary)
    if run['history']['linear_solves']:write_csv(REPORT/(name+'_linear.csv'),run['history']['linear_solves'])
    return result

def main():
    mode=argparse.ArgumentParser();mode.add_argument('mode',choices=['short','full']);args=mode.parse_args()
    assert load('petsc_smoke')['status']=='PASS','Smoke must pass before production setup'
    require_frozen_inputs()
    xml=ROOT/'configs/sv1_1/sv_flow_petsc.xml';old=ROOT/'configs/sv_flow.xml'
    if args.mode=='short':
        assert not xml.exists(),'Do not overwrite a frozen production candidate'
        tree=petsc_ls(ET.parse(old));ET.indent(tree);tree.write(xml,encoding='utf-8',xml_declaration=True)
        compare_production_xml(old,xml);validate_petsc_xml(xml)
        (REPORT/'xml_diff.txt').write_text(''.join(difflib.unified_diff(old.read_text().splitlines(True),xml.read_text().splitlines(True),fromfile='SV1 frozen',tofile='SV1.1 PETSc')))
        write_json(REPORT/'production_invariance.json',{'status':'PASS','semantic_differences_only_LS':True,'old_xml_sha256':sha256(old),'new_xml_sha256':sha256(xml),'options':OPTIONS})
    else:assert load('petsc_short_gate')['status']=='PASS','Short comparison must pass first'
    mesh=ROOT/'outputs/sv1_1/SV_MESH'
    if not mesh.exists():mesh.symlink_to(ROOT/'outputs/sv1/SV_MESH',target_is_directory=True)
    case=ROOT/'outputs/sv1_1'/('vascular_'+args.mode);case.mkdir(parents=True,exist_ok=True)
    assert not (case/'4-procs').exists(),'Initial-condition run requires an empty case directory'
    shutil.copyfile(xml,case/'solver.xml')
    compare_production_xml(old,case/'solver.xml')
    policy=json.loads((ROOT/'configs/time_policy.json').read_text());name='petsc_'+args.mode
    run=run_solver(case,name,policy['dt_s'],10 if args.mode=='short' else None)
    qc=measure(case,name,run)
    try:linear_gate(run);linear=True;error=None
    except ValueError as exc:linear=False;error=str(exc)
    try:nonlinear_gate(run['history']);nonlinear=True
    except ValueError:nonlinear=False
    last=qc['states'][-1] if qc['states'] else None
    base=load('flow_qc','sv1')
    mass_improves=bool(last and last['step']==10 and last['epsilon_mass']<base['epsilon_mass'])
    finite=bool(last and last['velocity_finite'] and last['pressure_finite'])
    complete=qc['completed_steps']==(10 if args.mode=='short' else 400)
    passed=linear and nonlinear and finite and complete and (mass_improves if args.mode=='short' else qc['steady_reached'])
    reason=None if passed else ('PETSC_LINEAR_SOLVE_FAIL' if not linear else 'NONLINEAR_CONVERGENCE_FAIL' if not nonlinear else 'SHORT_RUN_MASS_NOT_IMPROVED' if args.mode=='short' and not mass_improves else 'STEADY_NOT_REACHED')
    gate={'status':'PASS' if passed else 'FAIL','reason':reason,'linear_pass':linear,'linear_error':error,'nonlinear_pass':nonlinear,'finite_fields':finite,
          'complete':complete,'mass_improves':mass_improves if args.mode=='short' else None,'baseline_epsilon_mass':base['epsilon_mass'],'final_measurement':last}
    write_json(REPORT/(name+'_gate.json'),gate);print(json.dumps(gate,indent=2))
    return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
