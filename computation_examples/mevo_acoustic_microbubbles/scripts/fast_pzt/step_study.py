"""Only the six requested macro steps on one fixed formation scene."""
import argparse
import json
import numpy as np
import h5py
from .common import BASE,ROOT,write
from .experiment import run_case


def main(label,contact_iterations,candidates_us=None,reference=None):
    out=BASE/label; out.mkdir(parents=True,exist_ok=False)
    protocol=json.loads((ROOT/'configs/fast_pzt/f1_protocol.json').read_text()); write(out/'FROZEN_PROTOCOL.json',protocol)
    scene=protocol['formation']; rows=[]
    candidates=protocol['time_step_candidates_s'] if candidates_us is None else [v*1e-6 for v in candidates_us]
    for dt in candidates:
        cfg=dict(scene,dt_s=dt,events=[{'time_s':0,'pressure_Pa':scene['pressure_Pa']}],
                 physics={'contact_iterations':contact_iterations})
        name=f'dt_{round(dt*1e6)}us'; summary=run_case(out/name,cfg)
        rows.append(dict(name=name,**summary))
        write(out/'PARTIAL.json',dict(rows=rows))
    refrow=rows[0] if reference is None else json.loads((BASE/reference/'RESULT.json').read_text())['rows'][0]
    refdir=out if reference is None else BASE/reference
    with h5py.File(refdir/refrow['name']/'trajectory.h5') as h:
        ref={k:h['metrics/'+k][:] for k in ['cluster_count','largest_size','global_rg_m','mean_nearest_neighbor_m','cluster_sizes']}
    criteria=protocol['step_selection']
    for row in rows:
        with h5py.File(out/row['name']/'trajectory.h5') as h:
            values={k:h['metrics/'+k][:] for k in ref}
            if len(values['global_rg_m'])!=len(ref['global_rg_m']):
                row['step_statistics_pass']=False; row['comparison_reason']='Incomplete/common clocks unavailable'; continue
        errs={}
        for key in ['global_rg_m','mean_nearest_neighbor_m']:
            errs[key+'_relative_rms']=float(np.sqrt(np.mean((values[key]-ref[key])**2))/np.sqrt(np.mean(ref[key]**2)))
        errs['endpoint_cluster_count_difference']=int(abs(values['cluster_count'][-1,0]-ref['cluster_count'][-1,0]))
        errs['endpoint_largest_size_difference']=int(abs(values['largest_size'][-1,0]-ref['largest_size'][-1,0]))
        def distribution(sizes):
            sizes=sizes[sizes>0]; hist=np.bincount(sizes,minlength=scene['count']+1); return hist/hist.sum()
        errs['cluster_distribution_total_variation']=float(np.sum(abs(distribution(values['cluster_sizes'][-1,0])-distribution(ref['cluster_sizes'][-1,0])))/2)
        tref=refrow['formation_time_s'][0]; t=row['formation_time_s'][0]
        errs['formation_time_relative_difference']=None if t is None or tref is None else abs(t-tref)/tref
        row['comparison']=errs
        row['step_statistics_pass']=bool(
            errs['global_rg_m_relative_rms']<=criteria['global_rg_series_relative_rms'] and
            errs['mean_nearest_neighbor_m_relative_rms']<=criteria['mean_nearest_neighbor_series_relative_rms'] and
            errs['endpoint_cluster_count_difference']<=criteria['endpoint_cluster_count_absolute_difference'] and
            errs['endpoint_largest_size_difference']<=criteria['endpoint_largest_size_absolute_difference'] and
            errs['cluster_distribution_total_variation']<=criteria['endpoint_cluster_size_distribution_total_variation'] and
            errs['formation_time_relative_difference'] is not None and errs['formation_time_relative_difference']<=criteria['formation_time_relative_difference'])
    good=[r for r in rows if r['numerical_quality']=='PASS' and r.get('step_statistics_pass')]
    best=min(good,key=lambda r:r['wall_s']) if good else None
    write(out/'RESULT.json',dict(complete=True,rows=rows,selected=None if best is None else best['name'],
        selected_dt_s=None if best is None else best['nominal_dt_s'],contact_iterations=contact_iterations,
        reference_directory=str(refdir),reference_run=refrow['name']))


if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('--label',required=True); a.add_argument('--contact-iterations',type=int,default=4)
    a.add_argument('--candidates-us',type=float,nargs='+'); a.add_argument('--reference')
    q=a.parse_args(); main(q.label,q.contact_iterations,q.candidates_us,q.reference)
