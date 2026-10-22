"""Frozen paired experiments; same time, particles, F1/F2 and admissibility."""
import time
import numpy as np,h5py
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from .common import *
from .run import task,batch
from .slow import batch as slow_batch

def compare(t):
    a=read(BASE/'hf'/t['id']/'SUMMARY.json');b=read(BASE/'lf'/t['id']/'SUMMARY.json');row=dict(id=t['id'],HF_status=a['status'],LF_status=b['status'],state=t['state'],action=t['action'],duration_s=t['duration'],PASS=False)
    if a['status']!='PASS' or b['status']!='PASS':return row
    with h5py.File(BASE/a['trajectory_path']) as h,h5py.File(BASE/b['trajectory_path']) as z:
        ht=h['time_s'][:];xx=h['positions_m'];lt=z['sample_time_s'][:];yy=z['sample_positions_m'][:];ix=np.unique(np.rint(np.linspace(0,len(lt)-1,101)).astype(int));err=[];com=[];pair=[]
        for i in ix:
            j=int(np.argmin(abs(ht-lt[i])));dt=abs(ht[j]-lt[i]);assert dt<1e-12,'No interpolation of HF truth in validation'
            d=xx[j]-yy[i];err.append(np.linalg.norm(d,axis=1).max());com.append(np.linalg.norm(d.mean(0)));pair.append(np.max(abs(pdist(xx[j])-pdist(yy[i]))))
    af=a['final'];bf=b['final'];row.update(max_fixed_ID_m=max(err),max_COM_m=max(com),max_pair_distance_m=max(pair),Rg_error_m=abs(af['Rg_3d_m']-bf['Rg_3d_m']),structure_peak_error=abs(af['structure_peak_S']-bf['structure_peak_S']),stripe_descriptor_agrees=af['stripe_descriptor_pass']==bf['stripe_descriptor_pass'],lattice_descriptor_agrees=af['lattice_descriptor_pass']==bf['lattice_descriptor_pass'],HF_wall_s=a['task_wall_s'],LF_wall_s=b['wall_s'],HF_fallback_steps=b['HF_fallback_steps'],speedup=a['task_wall_s']/b['wall_s'],LF_min_whole_interval_margin_m=b['whole_interval_min_margin_m'],PASS=max(err)<=read(CONFIG)['LF']['fixed_ID_limit_m'])
    return row

def main():
    begin=time.perf_counter();assert read(BASE/'HISTORICAL_REPLAY.json')['bitwise_positions_equal']
    tasks=[task('pair_r0_'+a,a) for a in ['OFF','X_PAIR','Y_PAIR','XY_INPHASE','XY_QUADRATURE','ROT_PLUS','ROT_MINUS','X_SINGLE']]
    tasks += [task(f'pair_v{j}_{a}',a,state=f'N64_wide_v{j}') for j in range(3) for a in ['OFF','X_PAIR','ROT_PLUS','ROT_MINUS']]
    tasks += [task('pair_compact','X_PAIR','N64_compact_reference'),task('pair_contact','X_PAIR','N32_contact_transition',duration=.005)]
    tasks += [task(f'cost_N{n}','X_PAIR',f'N{n}_wide_r0',duration=.0002) for n in [32,128,200]]
    write(BASE/'PAIRED_PLAN.json',dict(tasks=tasks,identity=identity(),limits=read(CONFIG)['budget']))
    a=batch(tasks);b=slow_batch(tasks);rows=[compare(t) for t in tasks]
    pairs=[r for r in rows if 'max_fixed_ID_m' in r];rank=[]
    for key in ['structure_peak_S','relative_fixed_ID_rms_m']:
        u=[x['final'][key] for x,y in zip(a,b) if x['status']==y['status']=='PASS'];v=[y['final'][key] for x,y in zip(a,b) if x['status']==y['status']=='PASS'];rank.append(dict(metric=key,rho=float(spearmanr(u,v).statistic)))
    successes=sum(r['final']['stripe_descriptor_pass'] or r['final']['lattice_descriptor_pass'] for r in a if r['status']=='PASS')
    write(BASE/'PAIRED_VALIDATION.json',dict(rows=rows,ranks=rank,numeric_pass=all(r['PASS'] for r in rows),HF_pattern_successes=successes,scientific_screening='INCONCLUSIVE' if not successes else 'REQUIRES_NULL_AND_INDEPENDENT_CONFIRMATION',false_negatives=sum(x['status']=='PASS' and x['final']['stripe_descriptor_pass'] and (y['status']!='PASS' or not y['final']['stripe_descriptor_pass']) for x,y in zip(a,b)),wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
