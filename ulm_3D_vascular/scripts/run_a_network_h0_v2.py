"""Compute H0, sensitivities and pressure transfer; v1 and production remain read-only."""
from pathlib import Path
import sys, time, resource, json, platform
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.idealized_h0 import load_h0, solve_operating_point, downstream_resistance, external_radius_variant, MODEL_NAME, TERMINAL_MODEL
from network_1d0d.extension_transfer import measure_extensions, pressure_transfer
from network_1d0d.audit import write_json, write_csv, sha256
from network_1d0d.hydraulic_resistance import MU_PA_S, RHO_KG_M3, ROI_TARGET_Q_M3_S


def run(output):
    started=time.perf_counter();v1=ROOT/'reports/a_network_1d0d_boundary_v1';data=output/'data';data.mkdir(parents=True,exist_ok=True)
    domain=load_h0(v1);baseline=solve_operating_point(domain)
    write_json(data/'roi_port_flow_sign_convention.json',dict(edge_algebraic_direction_is_not_blood_flow=True,ports=domain.signs))
    np.savez_compressed(data/'analysis_A_H0_graph_si.npz',ids=domain.ids,xyz_m=domain.xyz_m,radius_m=domain.radius_m,edges=domain.edges,
                        original_edge=domain.raw_edge,fractions=domain.fractions,roi_internal_edge_mask=domain.internal_roi,
                        terminal_indices=domain.terminals,source_index=domain.source,port_indices=domain.port_nodes)
    write_csv(data/'unit_node_pressure.csv',[dict(node_id=int(n),pressure_Pa=float(p)) for n,p in zip(domain.ids,baseline['unit'].pressure)])
    write_csv(data/'unit_edge_flow.csv',[dict(edge_id=i,u=int(domain.ids[u]),v=int(domain.ids[v]),Q_u_to_v_m3s=float(baseline['unit'].flow[i])) for i,(u,v) in enumerate(domain.edges)])
    if baseline['status']!='PASS':
        write_json(data/'final_summary.json',dict(final_stage_status=baseline['status'],model_name=MODEL_NAME,unit_ROI_inlet_Q=float(baseline['unit_roi_q'][0])))
        return
    write_json(data/'network_mass_balance_audit.json',baseline['mass']|dict(solver=baseline['solver_audit']))
    write_json(data/'roi_operating_point_H0.json',dict(status=baseline['status'],ports=baseline['ports'],scaling_lambda=baseline['scaling_lambda'],unit_roi_signed_Q=baseline['unit_roi_q'],mass=baseline['mass']))
    np.savez_compressed(data/'H0_solution_si.npz',pressure_Pa=baseline['pressure'],edge_Q_m3s=baseline['flow'],node_outflow_m3s=baseline['node_outflow'],internal_node_indices=baseline['internal_nodes'])
    write_csv(data/'scaled_node_pressure.csv',[dict(node_id=int(n),pressure_Pa=float(p),signed_node_outflow_m3s=float(q)) for n,p,q in zip(domain.ids,baseline['pressure'],baseline['node_outflow'])])
    write_csv(data/'scaled_edge_flow.csv',[dict(edge_id=i,u=int(domain.ids[u]),v=int(domain.ids[v]),Q_u_to_v_m3s=float(baseline['flow'][i]),R_Pa_s_m3=float(domain.resistances()[i])) for i,(u,v) in enumerate(domain.edges)])
    write_csv(data/'terminal_outflows.csv',[dict(node_id=int(domain.ids[t]),pressure_Pa=float(baseline['pressure'][t]),outflow_m3s=float(-baseline['node_outflow'][t]),definition='IDEALIZED_REFERENCE_PRESSURE_TERMINAL') for t in domain.terminals])
    downstream={name:downstream_resistance(domain,baseline,name) for name in ['O1','O2','O3']}
    write_json(data/'downstream_resistance_audit.json',downstream)
    loo=[]
    for t in domain.terminals:
        try:
            variant=solve_operating_point(domain,removed_terminal=int(t))
            if variant['status']!='PASS':raise ValueError(variant['status'])
            row=dict(removed_terminal_id=int(domain.ids[t]),status='PASS',reason='',scaling_lambda=variant['scaling_lambda'],source_pressure_Pa=float(variant['pressure'][domain.source]),mass_residual=variant['mass']['global_relative_residual'])
            for i,name in enumerate(['O1','O2','O3']):row[name+'_fraction']=float(variant['fractions'][i]);row[name+'_difference_pp']=float(100*(variant['fractions'][i]-baseline['fractions'][i]))
        except (ValueError,RuntimeError) as exc:
            row=dict(removed_terminal_id=int(domain.ids[t]),status='INVALID_VARIANT',reason=str(exc),scaling_lambda=None,source_pressure_Pa=None,mass_residual=None)
            for name in ['O1','O2','O3']:row[name+'_fraction']=None;row[name+'_difference_pp']=None
        loo.append(row)
    write_csv(data/'terminal_leave_one_out_sensitivity.csv',loo)
    sensitivity={}
    valid=[row for row in loo if row['status']=='PASS']
    for i,name in enumerate(['O1','O2','O3']):
        values=np.array([row[name+'_fraction'] for row in valid]);delta=100*(values-baseline['fractions'][i]);worst=int(np.argmax(abs(delta)))
        sensitivity[name]=dict(baseline=float(baseline['fractions'][i]),min=float(values.min()),max=float(values.max()),P05_P50_P95=np.quantile(values,[.05,.5,.95]),
                               maximum_absolute_difference_pp=float(abs(delta[worst])),signed_worst_difference_pp=float(delta[worst]),most_influential_terminal_id=valid[worst]['removed_terminal_id'])
    write_json(data/'terminal_sensitivity_summary.json',dict(valid_count=len(valid),invalid_count=len(loo)-len(valid),ports=sensitivity,scope='One sink removed at a time, all other H0 BCs unchanged; rematch ROI Q each time'))
    radii_rows=[]
    for label,name,factor in [('global_'+str(f),'GLOBAL',f) for f in [.9,.95,1.,1.05,1.1]]+[(name+'_'+str(f),name,f) for name in ['O1','O2'] for f in [.95,1.05]]:
        radii=domain.radius_m*factor if name=='GLOBAL' else external_radius_variant(domain,name,factor)[0]
        variant=solve_operating_point(domain,radii=radii)
        row=dict(variant=label,scope=name,radius_factor=factor,status=variant['status'],source_pressure_Pa=float(variant['pressure'][domain.source]),source_pressure_ratio=float(variant['scaling_lambda']/baseline['scaling_lambda']),
                 ROI_inlet_pressure_Pa=variant['ports'][0]['pressure_realcut_Pa'],ROI_inlet_Q=variant['ports'][0]['signed_Q_m3s'],
                 perturbation_rule='all radii scaled' if name=='GLOBAL' else 'all saved exterior component node radii scaled except the real-cut node, which is held fixed to preserve ROI geometry')
        for i,p in enumerate(variant['ports'][1:]):row[p['port']+'_fraction']=p['signed_fraction'];row[p['port']+'_pressure_Pa']=p['pressure_realcut_Pa'];row[p['port']+'_difference_pp']=100*(p['signed_fraction']-baseline['fractions'][i])
        radii_rows.append(row)
    write_csv(data/'radius_hydraulic_sensitivity.csv',radii_rows)
    case=ROOT.parent/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps'
    extensions=measure_extensions(case,domain.ports,40)
    write_json(data/'current_FEM_extension_sections.json',extensions)
    transfer=pressure_transfer(baseline['ports'],extensions,sha256(v1/'data/roi_ports_in_a.json'),sha256(data/'H0_solution_si.npz'))
    write_json(data/'network_to_fem_pressure_transfer.json',transfer);write_json(data/'roi_fixed_pressure_bc_H0.json',transfer)
    current=np.array([.04257917,.85205026,.10537057]);frac=baseline['fractions']
    summary=dict(model_name=MODEL_NAME,model_version=2,hydraulic_source_node=2410,source_definition='IDEALIZED_STRUCTURAL_ROOT',terminal_model=TERMINAL_MODEL,
                 terminal_count=len(domain.terminals),p_ref_Pa=0.,mu_Pa_s=MU_PA_S,rho_kg_m3=RHO_KG_M3,ROI_target_Q_m3s=ROI_TARGET_Q_M3_S,
                 original_analysis_nodes=7419,original_analysis_edges=7418,solver_nodes=len(domain.ids),solver_edges=len(domain.edges),excluded_reference_components=42,
                 unit_source_pressure_Pa=1.,unit_ROI_inlet_Q=float(baseline['unit_roi_q'][0]),scaling_lambda=baseline['scaling_lambda'],scaled_source_pressure_Pa=float(baseline['pressure'][domain.source]),
                 A_total_source_inflow=baseline['mass']['A_total_source_inflow'],A_total_terminal_outflow=baseline['mass']['A_total_terminal_outflow'],network_mass_residual=baseline['mass'],
                 ROI_inlet_pressure_realcut=baseline['ports'][0]['pressure_realcut_Pa'],ROI_inlet_Q=baseline['ports'][0]['signed_Q_m3s'],current_3D_split=current,H0_split=frac,
                 split_difference=dict(percentage_points=100*(frac-current),relative_change=(frac-current)/current),terminal_sensitivity_range=sensitivity,radius_sensitivity=radii_rows,
                 downstream_resistances=downstream,new_3D_case_status='NOT_CREATED_YET',new_3D_split=None,final_stage_status='A_NETWORK_H0_BASELINE_READY_3D_PENDING',
                 source_pressure_interpretation='Idealized H0 gauge pressure scale at matched ROI flow, not measured mouse arterial pressure',
                 workers=1,server_used_for_network=False,network_runtime_seconds=time.perf_counter()-started,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,host=platform.node())
    for p in baseline['ports'][1:]:
        name=p['port'];summary[name+'_pressure_realcut']=p['pressure_realcut_Pa'];summary[name+'_Q']=p['signed_Q_m3s'];summary[name+'_fraction']=p['signed_fraction']
    for p in transfer['ports']:
        name=p['port'];summary[name+'_extension_R']=p['R_extension_Pa_s_m3'];summary[name+'_pressure_cap_raw']=p['pressure_cap_raw_Pa'];summary[name+'_pressure_cap_shifted']=p['pressure_cap_shifted_Pa']
    gates=dict(H0_sparse_solve=baseline['status']=='PASS',ROI_inlet_direction=baseline['ports'][0]['signed_Q_m3s']>0,network_mass_balance=baseline['mass']['status']=='PASS',
               exactly_four_ports=len(domain.ports)==4,cap_transfer_available=transfer['status']=='PASS',fixed_BC_finite=all(np.isfinite(p['pressure_cap_shifted_Pa']) for p in transfer['ports']),
               pressure_differences_preserved=transfer['max_pressure_difference_preservation_error_Pa']<1e-9)
    summary['pre_3D_gates']=gates
    write_json(data/'final_summary.json',summary)
    print(json.dumps(dict(status=summary['final_stage_status'],terminal_count=len(domain.terminals),H0_split=frac.tolist(),source_pressure=summary['scaled_source_pressure_Pa'],cap_pressures=[p['pressure_cap_shifted_Pa'] for p in transfer['ports']],runtime_s=summary['network_runtime_seconds'],gates=gates),indent=2))


if __name__=='__main__':run(ROOT/'reports/a_network_1d0d_boundary_v2_idealized')
