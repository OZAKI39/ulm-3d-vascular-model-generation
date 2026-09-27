"""Independent actual-field/log validation of the ONE frozen-design CFD.

Read-only numerical validation precedes frozen-flow export. No solver launch,
parameter update, pressure tuning, or change to legacy production results.
"""
from pathlib import Path
import argparse
import json
import math
import shutil
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.audit import sha256,write_json,write_csv
from network_1d0d.balance_design import load_frozen_design,metrics
from network_1d0d.fem_h0_case import audit_xml_change


def validate(case,source,flow_root,report):
    case,source,flow_root,report=map(Path,(case,source,flow_root,report))
    sys.path[:0]=[str(flow_root/'src'),str(flow_root/'scripts/sv13q')]
    import pyvista as pv
    from sv_validation.postprocess import SolutionMeasurements
    from sv_validation.sv13 import SteadyStopMonitor,stop_gate
    from sv_validation.sv13n import checkpoint_one_rank
    from sv_validation.sv11 import linear_gate,nonlinear_gate
    from flow_parser import parse_solver_log
    pre=json.loads((case/'final_preflight.json').read_text())
    frozen=load_frozen_design(case/'frozen_balance_design.yaml')
    assert sha256(case/'frozen_balance_design.yaml')==pre['frozen_design_sha256']
    assert sha256(case/'best_feasible_fem_handoff.json')==pre['handoff_sha256']
    fields=audit_xml_change((source/'run/solver.xml').read_bytes(),(case/'run/solver.xml').read_bytes(),
        [pre['applied_cap_pressure_pa'][p] for p in ('O1','O2','O3')])
    for name,digest in pre['configuration_diff']['unchanged_input_hashes'].items():
        assert sha256(case/name)==sha256(source/name)==digest,name
    for name,digest in pre['input_hashes'].items():assert sha256(case/name)==digest,name
    execution=json.loads((case/'reports/execution.json').read_text())
    dispatch=json.loads((case/'reports/one_scientific_case_dispatch.json').read_text())
    completion=json.loads((case/'reports/dispatch_completion.json').read_text())
    assert dispatch['dispatch_count']==completion['new_scientific_CFD_run_count']==1
    assert completion['runner_exit_code']==0 and execution['status']=='PASS' and execution['exit_code']==0
    assert execution['solver_sha256']==pre['solver_sha256'] and execution['PETSc_library_sha256']==pre['petsc_library_sha256']
    assert sha256(case/'run/solver.log')==execution['log_sha256']
    policy=json.loads((case/'policy.json').read_text())
    history=parse_solver_log((case/'run/solver.log').read_text(),policy['dt_s'])
    assert history==execution['history']
    linear_gate(execution);nonlinear_gate(history)
    measure=SolutionMeasurements(case/'SV_MESH/mesh_arrays.npz',policy['Q_target_m3_s'],policy['Umean_m_s'])
    monitor=SteadyStopMonitor(measure,policy);qualified=[];snapshots=[]
    for path in sorted((case/'run/1-procs').glob('result_*.vtu'),key=lambda p:int(p.stem.rsplit('_',1)[1])):
        step=int(path.stem.rsplit('_',1)[1])
        if step==0 or step%policy['save_interval_steps']:continue
        u,p=measure.read(path)
        m=dict(measure.measure(u,p),step=step,time_s=step*policy['dt_s'])
        if monitor.observe(m,u):qualified.append(step)
        snapshots.append(dict(step=step,path=str(path.relative_to(case)),sha256=sha256(path)))
    assert execution['stop']['step'] in qualified,'Independent steady sequence did not qualify'
    final=case/execution['final_vtu'];u,p=measure.read(final);grid=pv.read(final)
    assert np.array_equal(grid.points,measure.points)
    assert np.all(grid.celltypes==10)
    np.testing.assert_array_equal(np.sort(grid.cells.reshape(-1,5)[:,1:],axis=1),np.sort(measure.tetra,axis=1))
    m=measure.measure(u,p)
    assert stop_gate(monitor.intervals,m,0,0,policy)
    assert m['Q_in_m3_s']>0 and m['epsilon_Q']<=policy['mass_limit'] and m['epsilon_mass']<=policy['mass_limit']
    assert m['wall_noslip_pass'] and m['velocity_finite'] and m['pressure_finite']
    qualifying=final.parent/('result_%03d.vtu'%execution['stop']['step'])
    uq,pq=measure.read(qualifying);mq=measure.measure(uq,pq)
    eu=measure.velocity_l2(u-uq)/measure.velocity_l2(u)
    eq=max(abs(m['signed_outward_boundary_flows_m3_s'][n]-mq['signed_outward_boundary_flows_m3_s'][n])/measure.Q for n in m['signed_outward_boundary_flows_m3_s'])
    assert eu<=policy['velocity_change_limit'] and eq<=policy['flow_change_limit']
    checkpoint=final.with_name('stFile_%03d.bin'%execution['final_step'])
    cp=checkpoint_one_rank(checkpoint,execution['final_step'],policy['dt_s'])
    roles=('OUTLET_01','OUTLET_02','OUTLET_03');ports=('O1','O2','O3')
    f=np.array([m['outlet_flows_m3_s'][n]/m['Q_in_m3_s'] for n in roles])
    assert np.all(f>=0),'Outlet backflow found; freeze blocked'
    m['signed_fractions_Qout_over_Qin']=dict(zip(ports,f.tolist()))
    m['port_backflow']={p:bool(v<0) for p,v in zip(ports,f)}
    m['epsilon_mass_compensated']=abs(math.fsum([*m['outlet_flows_m3_s'].values(),-m['Q_in_m3_s']]))/measure.Q
    old_path=source/'frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu'
    old_manifest=json.loads((source/'frozen_flow/manifest.json').read_text())
    assert sha256(old_path)==old_manifest['files'][old_path.name]['sha256']
    old_u,old_p=measure.read(old_path);old=measure.measure(old_u,old_p)
    old_f=np.array([old['outlet_flows_m3_s'][n]/old['Q_in_m3_s'] for n in roles])
    zero_f=np.array([frozen['optimized']['prediction']['outlet_flow_fraction'][n] for n in ports])
    fm,om=metrics(f),metrics(old_f)
    comparison=dict(status='BALANCE_IMPROVED_IN_3D' if fm['J_balance']<om['J_balance'] and fm['range']<om['range'] else 'BALANCE_NOT_IMPROVED_IN_3D',
        old_H0_3D_fractions=old_f,new_3D_fractions=f,design_0D_fractions=zero_f,
        discrepancy_3D_minus_0D=f-zero_f,discrepancy_percentage_points=100*(f-zero_f),
        old_H0_3D_metrics=om,new_3D_metrics=fm,
        delta_J=om['J_balance']-fm['J_balance'],relative_J_reduction=1-fm['J_balance']/om['J_balance'],
        range_reduction=om['range']-fm['range'],std_reduction=om['std']-fm['std'],
        old_native_flow_sha256=sha256(old_path),CFD_based_retuning=False,
        interpretation='Independent 3D forward result; equal split or 0D agreement are not acceptance gates')
    out=case/'frozen_flow';out.mkdir(exist_ok=False)
    exported=out/'steady_flow_mean_2p0_mmps_A_best_feasible_balance.vtu'
    shutil.copyfile(final,exported);shutil.copyfile(checkpoint,out/checkpoint.name)
    np.savez_compressed(out/'flow_arrays_si.npz',points_m=measure.points,tetra=measure.tetra,
        boundary_triangles=measure.boundary,facet_tags=measure.tags,velocity_m_s=u,pressure_pa=p)
    manifest=dict(status='PASS',role='Independent one-way CFD prediction from FROZEN_DESIGN',
        particle_production_promoted=False,new_scientific_CFD_run_count=1,CFD_based_retuning=False,
        design_envelope=frozen['design_bounds'],optimized_0D_parameters=frozen['optimized']['parameters'],
        design_0D_prediction=frozen['optimized']['prediction'],design_0D_balance_metrics=frozen['optimized']['metrics'],
        FEM_cap_pressures_pa=pre['applied_cap_pressure_pa'],actual_3D_fractions=dict(zip(ports,f.tolist())),
        actual_3D_balance_metrics=fm,discrepancy_3D_minus_0D=dict(zip(ports,(f-zero_f).tolist())),
        old_H0_comparison=comparison,source_frozen_design_sha256=pre['frozen_design_sha256'],
        geometry_sha256=sha256(case/'SV_MESH/mesh-complete.exterior.vtp'),mesh_sha256=sha256(case/'SV_MESH/mesh-complete.mesh.vtu'),
        wall_sha256=sha256(case/'SV_MESH/mesh-surfaces/WALL.vtp'),solver_XML_sha256=sha256(case/'run/solver.xml'),
        solver_binary_sha256=execution['solver_sha256'],native_output_sha256=sha256(final),
        final_step=execution['final_step'],node_count=len(measure.points),tetra_count=len(measure.tetra),
        units=dict(points='m',velocity='m/s',pressure='Pa',flow='m^3/s',WSS='Pa'),
        files={p.name:dict(sha256=sha256(p),bytes=p.stat().st_size) for p in out.iterdir()})
    write_json(out/'manifest.json',manifest)
    result=dict(status='PASS',measurements=m,comparison=comparison,configuration_diff=fields,
        independent_log_history_matches=True,linear_gate='PASS',nonlinear_gate='PASS',
        independently_recomputed_intervals=monitor.intervals,independently_qualifying_steps=qualified,
        independent_snapshots=snapshots,final_to_stop_relative_velocity_change=eu,final_to_stop_normalized_flow_change=eq,
        checkpoint=cp,exit_code=execution['exit_code'],wall_time_s=execution['wall_time_s'],final_step=execution['final_step'],
        frozen_flow=str(exported.resolve()),frozen_flow_sha256=sha256(exported),
        new_scientific_CFD_run_count=1,CFD_based_retuning=False,particle_RBC_calls=0,
        internal_section_or_mesh_independence_claim=False)
    write_json(case/'reports/independent_validation.json',result)
    write_json(report/'independent_validation.json',result)
    write_csv(report/'balance_comparison.csv',[dict(case=label,**{p:float(v) for p,v in zip(ports,frac)},**metrics(frac))
        for label,frac in [('old_H0_3D',old_f),('best_feasible_0D',zero_f),('best_feasible_3D',f)]])
    write_csv(report/'final_3D_ports.csv',[dict(port=p,flow_m3_s=m['outlet_flows_m3_s'][role],fraction=float(f[i]),
        area_average_pressure_pa=m['area_average_pressure_pa'][role],difference_from_0D=float(f[i]-zero_f[i])) for i,(p,role) in enumerate(zip(ports,roles))])
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ('case','source-case','flow-root','report'):parser.add_argument('--'+k,type=Path,required=True)
    a=parser.parse_args()
    try:
        r=validate(a.case,a.source_case,a.flow_root,a.report)
        print(json.dumps({k:r[k] for k in ('status','exit_code','wall_time_s','final_step','measurements','comparison','frozen_flow_sha256')},indent=2,default=lambda a:a.tolist()))
    except Exception as e:
        write_json(a.report/'independent_validation_failure.json',dict(status='FAIL',error=type(e).__name__+': '+str(e),CFD_retry_permitted=False))
        raise
