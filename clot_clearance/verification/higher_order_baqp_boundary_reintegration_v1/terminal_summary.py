"""Only the requested terminal summary; no numerical work."""
import csv,json
from pathlib import Path
import numpy as np
V=Path(__file__).resolve().parent
def js(n):return json.loads((V/n).read_text())
def rows(n):
    with (V/n).open() as f:return list(csv.DictReader(f))
def main():
    lines=[]
    def out(s):lines.append(s)
    cfg=js('PD2_BOUNDARY_CONFIGURATION.json');status=js('FINAL_STATUS.json');target=js('PD2_TARGET_SUMMARIES.json');res=js('RESOURCE_SUMMARY.json')
    out('branch: clot-clearance/pfpd-progress-20261005 / HEAD: '+(V/'GIT_HEAD.txt').read_text().strip())
    out(f"N19: h={cfg['h_m']} m; delta={cfg['delta_m']} m; delta/h={cfg['delta_m']/cfg['h_m']}")
    out('PD2 mechanics: '+status['PD2_VARIATIONAL_MECHANICS'])
    out('Force stencil: 3 graph expansions; D0/D1/D2/D3=1/257/1797/5821; N27; fully closed YES')
    for a in js('PD2_AFFINE_BULK_GATE.json'):out(f"Affine bulk {a['kind']}: {a['bulk_force_scaled']:.9e}")
    bulk=rows('PD2_BULK_ORACLE.csv')
    for k in cfg['solve_order'][:4]:
        a=js('PD2_'+k.upper()+'.json');t=target[k]
        out(k+':')
        if 'load_checks' in t:
            c=t['load_checks'];out(f"  resultant={c['resultant_error']:.9e}; moment={c['moment_error']:.9e}; work={c['virtual_work_error']:.9e}")
        else:
            q=np.array([[float(r['q_'+c+'_N_m3']) for c in 'xyz'] for r in bulk if r['case']==k])
            out(f"  max |q_bulk|={np.linalg.norm(q,axis=1).max():.9e} N/m3; max |e_model|={t['bulk_model_max_N_m3']:.9e} N/m3; nodewise closure={t['nodewise_decomposition_error']:.9e}")
        out(f"  solver residual={a['force_residual']:.9e}; max Fbar={a['fields']['Fbar']['F_error_max_scaled']:.9e}; max Ftilde={a['fields']['Ftilde']['F_error_max_scaled']:.9e}; max Psi+={max(r['Psi_plus_error_max_scaled'] for r in a['fields'].values()):.9e}")
    c=js('PD2_CLASSICAL_SHEAR.json');e=js('PD2_CLASSICAL_TARGET_RESIDUAL.json')
    out(f"Classical shear: target max |e_model V|={e['model_force_max_N']:.9e} N; identity error={e['error_scaled']:.9e}; final residual={c['force_residual']:.9e}")
    out(f"  displacement shift={c['displacement_shift_max_m']:.9e} m; Fbar shift={c['fields']['Fbar']['F_error_max']:.9e}; Ftilde shift={c['fields']['Ftilde']['F_error_max']:.9e}; Psi+ shift={max(r['Psi_plus_error_max_Pa'] for r in c['fields'].values()):.9e} Pa")
    comparison=rows('C1_PD2_BOUNDARY_COMPARISON.csv')
    for k in cfg['solve_order'][:4]:
        old=next(a for a in comparison if a['case']==k and a['backend']=='C1');new=next(a for a in comparison if a['case']==k and a['backend']=='PD2')
        out(f"Support {k}: d50/d90/d95/d99="+'/'.join(new[f'd{x}_h'] for x in [50,90,95,99])+' h')
        out(f"  magnitude sum PD2/C1={float(new['L1_N'])/float(old['L1_N']):.9e}; d99 C1/PD2={old['d99_h']}/{new['d99_h']} h; max load difference={float(new['C1_minus_PD2_max_N']):.9e} N")
    for a in js('PD2_CLOSED_CENTER_LOCALIZATION.json'):
        if a['fully_closed']:out(f"Closed-center {a['case']}: {a['boundary_force_scaled']:.9e}")
    for k in ['HIGHER_ORDER_BOUNDARY_MECHANICS','CURVED_CLOT_WALL_READINESS','PRODUCTION_BAQP','PHASE_FORCE','C5_C7']:out(k+': '+status[k])
    out('Targeted tests: 12 passed; full pytest NOT RUN')
    out(f"CPU workers: 8 oracle / up to 2 equilibrium after measured memory proof; BLAS=1; wall={res['scientific_wall_seconds']:.3f} s; peak RAM={res['peak_process_tree_RSS_bytes']/2**30:.3f} GiB; minimum headroom={res['minimum_available_RAM_fraction']:.2%}; GPU_NOT_USED_BY_DESIGN")
    out('REVIEW: '+str(V/'REVIEW_HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1.md'))
    out('Figures: '+str(V.parents[1]/'results/higher_order_baqp_boundary_reintegration_v1/figures'))
    text='\n'.join(lines)+'\n';(V/'TERMINAL_SUMMARY.txt').write_text(text);print(text)
if __name__=='__main__':main()
