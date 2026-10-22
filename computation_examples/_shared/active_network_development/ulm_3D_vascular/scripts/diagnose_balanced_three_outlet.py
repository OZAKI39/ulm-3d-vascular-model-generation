"""Diagnose a REJECTED equal-split fit; never freeze, create or launch CFD.

Consumes the one recorded fit. Only additional full-network 0D forwards are
used to check its limiting behavior. No optimizer rerun or parameter update.
"""
from pathlib import Path
from dataclasses import replace
from collections import Counter
import argparse
import csv
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from network_1d0d.audit import write_json,write_csv,sha256
from network_1d0d.parameter_config import read_configuration,code_provenance
from network_1d0d.parameterized_hydraulics import load_geometry_cache,solve_parameterized_operating_point
from network_1d0d.balanced_feasibility import balance_metrics,passive_o1_o3_certificate,require_balanced_fit


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--fit-directory',type=Path,required=True)
    parser.add_argument('--report-directory',type=Path,required=True)
    args=parser.parse_args();fit=args.fit_directory;report=args.report_directory
    loaded=read_configuration(args.config)
    rejection=json.loads((fit/'fit_rejection.json').read_text())
    if rejection['status']!='REJECTED_UNIDENTIFIABLE':raise ValueError('This diagnostic consumes an existing rejected fit only')
    if rejection['provenance']['code_sha256']!=code_provenance():raise ValueError('Fit model source changed')
    if rejection['provenance']['source_config_sha256']!=loaded['config_sha256']:raise ValueError('Fit config changed')
    trace=list(csv.DictReader((fit/'parameter_trace.csv').open()))
    first,last=trace[0],trace[-1]
    assert last['phase']=='final' and all(r['status']=='PASS' for r in trace)
    cache=load_geometry_cache(loaded['geometry_path'],loaded['ports_path'])
    cert=passive_o1_o3_certificate(cache,mu_pa_s=loaded['spec'].mu_pa_s,target_roi_flow_m3_s=loaded['target_q'])
    assert cert['status']=='EQUAL_SPLIT_IMPOSSIBLE_UNDER_PRESCRIBED_MODEL'
    spec=replace(loaded['spec'],s_O1=float(last['s_O1_dimensionless']),s_O2=float(last['s_O2_dimensionless']))
    state=solve_parameterized_operating_point(cache,spec,loaded['target_q'])
    predicted=state.summary()
    f=state.signed_outlet_fractions
    np.testing.assert_allclose(f,[float(last['f_'+p+'_dimensionless']) for p in ('O1','O2','O3')],rtol=1e-10,atol=1e-13)
    positive=bool(np.isfinite(state.effective_radius_m).all() and (state.effective_radius_m>0).all()
                  and np.isfinite(state.edge_resistance_pa_s_m3).all() and (state.edge_resistance_pa_s_m3>0).all())
    assert positive
    projection=np.array(cert['equal_sigma_least_squares_boundary_projection'])
    cert['last_fit_vs_projection_max_abs_fraction']=float(np.max(abs(f-projection)))
    assert cert['last_fit_vs_projection_max_abs_fraction']<1e-7
    summary=dict(status='BALANCED_0D_TARGET_NOT_REACHED',model_name=predicted['model_name'],
        parameter_source='SYNTHETIC_TARGET',purpose='SYNTHETIC_BALANCED_FLOW_DESIGN',
        parameters_are_accepted=False,initial_parameters=dict(s_O1=1.,s_O2=1.),
        final_parameters=predicted['parameters'],prediction=predicted,
        initial_residual=[float(first['flow_residual_'+p+'_dimensionless']) for p in ('O1','O2','O3')],
        final_residual=[float(last['flow_residual_'+p+'_dimensionless']) for p in ('O1','O2','O3')],
        initial_objective_norm=float(first['objective_norm_dimensionless']),
        final_objective_norm=float(last['objective_norm_dimensionless']),
        initial_sum_squared_residual=float(first['objective_norm_dimensionless'])**2,
        final_sum_squared_residual=float(last['objective_norm_dimensionless'])**2,
        nfev=None,nfev_note='Base fitter raised after final rank check before serializing SciPy status/nfev; not invented or inferred from total residual calls.',
        optimizer_status=dict(success=None,termination='UnidentifiableParameterizationError: final Jacobian rank deficient'),
        total_forward_evaluations=len(trace),evaluation_phase_counts=dict(Counter(r['phase'] for r in trace)),
        identifiability_audit=rejection['audit'],balance_metrics=balance_metrics(f),
        positive_finite_radius_and_resistance=positive,
        raw_rejection_label_caveat='The base fitter labels any final data-rank deficiency PRIOR_REGULARIZED_NOT_DATA_IDENTIFIED. Here prior_constraints=0: no prior was used. Interpret the explicit rejection/rank instead.',
        provenance=rejection['provenance'],input_hashes=dict(cache.input_hashes),
        source_geometry_sha256=cache.geometry_sha256,source_rejection_sha256=sha256(fit/'fit_rejection.json'),
        source_trace_sha256=sha256(fit/'parameter_trace.csv'),
        diagnostic_reconstruction='Final rejected trace parameters, one unchanged full-network forward; not another fit.',
        CFD_calls=0,particle_simulation_calls=0,CFD_feedback_used=False)
    try:require_balanced_fit(summary)
    except ValueError as exc:summary['freeze_and_CFD_gate']=dict(status='BLOCKED',reason=str(exc))
    else:raise AssertionError('Rejected design must never pass')
    write_json(fit/'fit_summary.json',summary)
    np.savez_compressed(fit/'network_state_si.npz',status='REJECTED_DIAGNOSTIC_ONLY_NOT_FROZEN',
        original_node_ids=cache.node_ids,pressure_pa=state.pressure_pa,edge_nodes=state.edge_nodes,
        edge_flow_m3_s=state.edge_flow_m3_s,node_outflow_m3_s=state.node_outflow_m3_s,
        edge_resistance_pa_s_m3=state.edge_resistance_pa_s_m3,
        original_node_radius_m=state.effective_radius_m,port_indices=cache.port_indices,
        port_names=np.array(['INLET','O1','O2','O3']))
    write_json(report/'feasibility_certificate.json',cert)
    rows=[]
    for scale in (1.,2.,5.,10.,100.,10000.):
        result=solve_parameterized_operating_point(cache,replace(spec,s_O1=scale),loaded['target_q'])
        r1=cert['paths']['O1']['resistance_pa_s_m3'];r3=cert['paths']['O3']['resistance_pa_s_m3']
        identity=result.port_pressure_pa[1]-result.port_pressure_pa[3]-(r3*result.port_flow_m3_s[3]-r1*result.port_flow_m3_s[1])
        assert abs(identity)<1e-7
        assert result.signed_outlet_fractions[0]<=cert['ratio_bound_Q_O1_over_Q_O3']*result.signed_outlet_fractions[2]+1e-11
        rows.append(dict(role='FULL_A_0D_DIAGNOSTIC_NO_OPTIMIZATION',s_O1_dimensionless=scale,s_O2_dimensionless=spec.s_O2,
            **{p+'_fraction_dimensionless':float(x) for p,x in zip(('O1','O2','O3'),result.signed_outlet_fractions)},
            P_O1_pa=float(result.port_pressure_pa[1]),fixed_path_identity_error_pa=float(identity),
            mass_residual_dimensionless=result.mass_audit['max_relative_residual']))
    write_csv(report/'saturation_checks.csv',rows)
    write_csv(report/'last_rejected_ports.csv',[dict(case='REJECTED_BALANCED_DESIGN',port=p,
        pressure_pa=float(pressure),flow_m3_s=float(q),fraction_of_inlet=float(q/state.port_flow_m3_s[0]))
        for p,pressure,q in zip(('INLET','O1','O2','O3'),state.port_pressure_pa,state.port_flow_m3_s)])
    progress=dict(status='BALANCED_0D_TARGET_NOT_REACHED',reason='Fixed ROI geometry and passive O1 downstream forbid equal O1/O3 flows with legacy-reference O3.',
        base_commit='f43c09b1e55cbd9702a65c600a464154fcf346f3',fit_workflows_run=1,
        fitting_forward_evaluations=len(trace),additional_diagnostic_0D_forwards=1+len(rows),
        stages={'1_base_commit':'PASS','2_old_tests':'PASS_103','3_balanced_config':'PASS',
          '4_pure_0D_fit':'EXECUTED_REJECTED','5_balanced_acceptance':'FAIL_TARGET_UNREACHABLE',
          **{s:'BLOCKED_BY_REQUIRED_0D_ACCEPTANCE' for s in ('6_freeze','7_handoff','8_create_case','9_CFD_preflight','10_CFD','11_3D_validation','12_freeze_flow','13_3D_comparison')}},
        scientific_CFD_cases_created=0,new_scientific_CFD_run_count=0,
        CFD_based_parameter_retuning=False,particle_RBC_run_count=0,
        server_state_changes=False,automatic_CFD_launch_registered=False,
        no_freeze_no_handoff_no_CFD=True)
    write_json(report/'workflow_status.json',progress)
    write_json(report/'diagnostic_provenance.json',dict(script_sha256=sha256(__file__),
        certificate_module_sha256=sha256(ROOT/'network_1d0d/balanced_feasibility.py'),
        input_config_sha256=loaded['config_sha256'],original_fit_rejection_sha256=sha256(fit/'fit_rejection.json'),
        code_sha256=code_provenance(),protected_core_code_unchanged_since_fit=True))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.3,6.3),layout='constrained')
    x=np.linspace(0,1,601);k=cert['ratio_bound_Q_O1_over_Q_O3']
    ax.fill_between(x,0,np.minimum(k*x,1-x),color='#c8dfef',label='Necessary passive-network region (superset)')
    ax.plot(x,np.minimum(k*x,1-x),color='#216b9b',lw=1.5)
    ax.plot([0,1],[1,0],':',color='gray',label='f(O1) + f(O3) = 1')
    ax.scatter([1/3],[1/3],s=130,marker='x',color='#c13f34',label='Requested equal-split target',zorder=5)
    ax.scatter([f[2]],[f[0]],s=70,color='#d28c1c',label='Last rejected fit (not frozen)',zorder=6)
    ax.scatter([float(first['f_O3_dimensionless'])],[float(first['f_O1_dimensionless'])],s=60,color='#253858',label='Initial H0 state',zorder=6)
    ax.annotate(f'f(O1) <= {k:.6f} f(O3)',(.65,k*.65),xytext=(.46,.4),arrowprops={'arrowstyle':'->'},fontsize=11)
    ax.set(xlim=(0,1),ylim=(0,.65),xlabel='O3 flow fraction (dimensionless)',ylabel='O1 flow fraction (dimensionless)',
        title='Balanced design is infeasible under the prescribed 0D model\nFixed ROI geometry; passive O1 downstream; O3 pressure = 0 Pa')
    ax.legend(loc='upper right',fontsize=9);ax.grid(alpha=.2)
    fig.savefig(report/'balanced_target_feasibility.png',dpi=220);fig.savefig(report/'balanced_target_feasibility.pdf');plt.close(fig)
    print(json.dumps(dict(status=progress['status'],last_candidate=predicted['parameters'],fractions=f.tolist(),
        balance=summary['balance_metrics'],ratio_bound=k,equal_split_required_P_O1_pa=cert['required_P_O1_minus_Pd_for_equal_pa'],
        CFD_runs=0),indent=2))


if __name__=='__main__':main()
