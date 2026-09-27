"""Independently validate the actual H0 FEM solution, without changing production."""
from pathlib import Path
import argparse
import sys
import json
import math
import shutil
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from network_1d0d.audit import sha256,write_json
from network_1d0d.fem_h0_case import audit_xml_change


def run(case,flow_root,report):
    sys.path[:0]=[str(flow_root/'src'),str(flow_root/'scripts/sv13q')]
    import pyvista as pv
    from sv_validation.postprocess import SolutionMeasurements
    from sv_validation.sv13 import stop_gate
    from sv_validation.sv13n import checkpoint_one_rank
    from sv_validation.sv11 import linear_gate,nonlinear_gate
    from flow_parser import parse_solver_log
    data=report/'data';old=flow_root/'flow_cases/mean-2p0-mmps'
    bc=json.loads((data/'roi_fixed_pressure_bc_H0.json').read_text())
    values=[r['pressure_cap_shifted_Pa'] for r in bc['ports']]
    diff=json.loads((data/'3d_case_config_diff.json').read_text())
    assert sha256(old/'run/solver.xml')==diff['original_solver_sha256']
    assert sha256(data/'roi_fixed_pressure_bc_H0.json')==diff['source_BC_SHA']
    fields=audit_xml_change((old/'run/solver.xml').read_bytes(),(case/'run/solver.xml').read_bytes(),values)
    for name,digest in diff['unchanged_input_hashes'].items():
        assert sha256(old/name)==sha256(case/name)==digest, name
    for name,digest in json.loads((case/'input_hashes.json').read_text()).items():
        assert sha256(case/name)==digest, name
    policy=json.loads((case/'policy.json').read_text())
    execution=json.loads((case/'reports/execution.json').read_text())
    assert execution['status']=='PASS' and execution['exit_code']==0
    assert sha256(case/'run/solver.log')==execution['log_sha256']
    # Parse the actual copied native log again, independently of stored metrics.
    history=parse_solver_log((case/'run/solver.log').read_text(),policy['dt_s'])
    assert history==execution['history']
    linear_gate(execution);nonlinear_gate(history)
    measure=SolutionMeasurements(case/'SV_MESH/mesh_arrays.npz',policy['Q_target_m3_s'],policy['Umean_m_s'])
    final=case/execution['final_vtu'];u,p=measure.read(final);grid=pv.read(final)
    assert np.array_equal(grid.points,measure.points)
    cells=grid.cells.reshape(-1,5)
    assert np.all(cells[:,0]==4) and np.all(grid.celltypes==10)
    assert np.array_equal(np.sort(cells[:,1:],axis=1),np.sort(measure.tetra,axis=1))
    m=measure.measure(u,p)
    compensated=math.fsum([*m['outlet_flows_m3_s'].values(),-m['Q_in_m3_s']])
    m.update(global_signed_flux_residual_compensated_m3s=compensated,
             epsilon_mass_compensated_boundary_sum=abs(compensated)/measure.Q,
             mass_roundoff_note='Original production gate uses ordinary double sums; compensated sum of the same signed boundary fluxes is also retained, with no forced zeros.')
    assert m['Q_in_m3_s']>0 and m['epsilon_Q']<=1e-6 and m['epsilon_mass']<=1e-6
    assert m['velocity_finite'] and m['pressure_finite'] and m['wall_noslip_pass']
    assert stop_gate(execution['intervals'],m,0,0,policy)
    qualified=execution['stop']['step']
    prior=next(x for x in final.parent.glob('result_*.vtu') if int(x.stem.rsplit('_',1)[1])==qualified)
    uq,pq=measure.read(prior);mq=measure.measure(uq,pq)
    final_u_change=measure.velocity_l2(u-uq)/measure.velocity_l2(u)
    final_q_change=max(abs(m['signed_outward_boundary_flows_m3_s'][r]-mq['signed_outward_boundary_flows_m3_s'][r])/measure.Q for r in m['signed_outward_boundary_flows_m3_s'])
    assert final_u_change<=policy['velocity_change_limit'] and final_q_change<=policy['flow_change_limit']
    checkpoint=final.with_name('stFile_%03d.bin'%execution['final_step'])
    cp=checkpoint_one_rank(checkpoint,execution['final_step'],policy['dt_s'])
    uold,pold=measure.read(old/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    old_m=measure.measure(uold,pold)
    names=['OUTLET_01','OUTLET_02','OUTLET_03']
    # Fractions explicitly use actual Qin; retain signed results even with backflow.
    split=np.array([m['outlet_flows_m3_s'][n]/m['Q_in_m3_s'] for n in names])
    current=np.array([old_m['outlet_flows_m3_s'][n]/old_m['Q_in_m3_s'] for n in names])
    summary=json.loads((data/'final_summary.json').read_text());h0=np.array(summary['H0_split'])
    positive=np.maximum(split,0);positive/=positive.sum()
    m.update(signed_fractions_Qout_over_Qin=split,positive_outflow_normalized_fractions=positive,
             port_backflow={n:m['outlet_flows_m3_s'][n]<0 for n in names},
             inlet_actual_mean_m_s=m['Q_in_m3_s']/policy['A_in_m2'])
    comparison=dict(current_3D_remeasured_split=current,provided_current_3D_split=summary['current_3D_split'],
                    H0_split=h0,new_3D_split=split,new_minus_H0_pp=100*(split-h0),
                    new_minus_current_3D_pp=100*(split-current),
                    new_vs_current_relative_change=(split-current)/current,
                    split_denominator='signed Qout / measured Qin; Qout/sum(Qout) also retained in measurements',
                    agreement_is_not_acceptance_gate=True)
    frozen=case/'frozen_flow';frozen.mkdir(exist_ok=True)
    export=frozen/'steady_flow_mean_2p0_mmps_A_H0.vtu'
    if export.exists(): assert sha256(export)==sha256(final), 'No conflicting frozen export overwrite'
    else: shutil.copyfile(final,export)
    shutil.copyfile(checkpoint,frozen/checkpoint.name)
    np.savez_compressed(frozen/'flow_arrays_si.npz',points_m=measure.points,tetra=measure.tetra,
                        boundary_triangles=measure.boundary,facet_tags=measure.tags,velocity_m_s=u,pressure_pa=p)
    assert sha256(export)==sha256(final)
    manifest=dict(case=case.name,status='PASS',role='Independent CFD result for review',particle_production_promoted=False,
                  formulation='P1/P1 + VMS unchanged',mesh_sha256=sha256(case/'SV_MESH/mesh-complete.mesh.vtu'),
                  geometry_sha256=sha256(case/'SV_MESH/mesh-complete.exterior.vtp'),
                  wall_sha256=sha256(case/'SV_MESH/mesh-surfaces/WALL.vtp'),
                  input_configuration_sha256=sha256(case/'run/solver.xml'),
                  raw_solver_vtu=str(final),final_step=execution['final_step'],node_count=len(measure.points),tetra_count=len(measure.tetra),
                  files={f.name:dict(sha256=sha256(f),bytes=f.stat().st_size) for f in frozen.iterdir() if f.is_file() and f.name!='manifest.json'},
                  units=dict(points='m',pressure='Pa',velocity='m/s',flow='m3/s'))
    write_json(frozen/'manifest.json',manifest)
    result=dict(status='PASS',measurements=m,comparison=comparison,policy=policy,configuration_diff=fields,
                unchanged_geometry_mesh_wall=True,native_linear_and_nonlinear_convergence='PASS',
                steady=dict(first_qualifying_step=qualified,final_step=execution['final_step'],
                            final_relative_velocity_change=final_u_change,final_normalized_flow_change=final_q_change,
                            final_five_intervals=execution['intervals'][-5:]),
                pressure_BC_type='Neumann pressure traction, not pointwise Dirichlet pressure; area mean p may differ by viscous normal traction',
                native_checkpoint=cp,export=manifest,wall_time_s=execution['wall_time_s'],
                MPI_ranks=execution['MPI_ranks'],OMP_NUM_THREADS=execution['OMP_NUM_THREADS'],
                peak_tree_RSS_MiB_sampled=execution['peak_tree_RSS_MiB_sampled'],
                peak_single_process_HWM_MiB=execution['peak_single_process_HWM_MiB'],
                child_rusage_maxrss_MiB=execution['child_rusage_maxrss_MiB'],
                internal_arbitrary_section_gate_applied=False)
    write_json(case/'reports/physics_validation_H0.json',result)
    write_json(data/'3d_physics_validation.json',result)
    summary.update(new_3D_case_status='PASS',new_3D_split=split,final_stage_status='A_NETWORK_H0_3D_VALIDATED',
                   new_3D_comparison=comparison,new_3D_mass_residual=m['epsilon_mass'],
                   new_3D_runtime_seconds=execution['wall_time_s'],new_3D_peak_RSS_MiB=execution['peak_tree_RSS_MiB_sampled'],
                   new_3D_workers=1,new_3D_server_used=True)
    write_json(data/'final_summary.json',summary)
    print(json.dumps(dict(status=result['status'],measurements=m,comparison=comparison),indent=2,default=lambda x:x.tolist()))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--flow-root',type=Path,required=True)
    parser.add_argument('--report',type=Path,default=ROOT/'reports/a_network_1d0d_boundary_v2_idealized')
    args=parser.parse_args();run(args.case,args.flow_root,args.report)
