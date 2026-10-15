"""Bounded production campaign, gated by both frozen strength calibrations."""
from pathlib import Path
import csv,json,time,math,hashlib
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .config import dump,PATTERNS
from .coherent_failure_card import OUTPUT,RESULTS
from .interface_model import prepare_interface_geometry,write_typed_data
from .current_position_provider import prepare_provider,lammps_provider,spatial_shape,pulse
from .calibrate_failure import prelude,dumps
from .runner import run_lammps
from .fragment_analysis import read_dump,_events
from .fragment_semantics import component_metrics, classify_detached_v2
from .failure_classifier import classify,clot_components
from .bond_utilization import utilization

def csv_out(path,rows):
    if not rows:path.write_text('case_id,status\n');return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def case_request(config,pattern,amplitude,duration_factor,stage,index):
    T=round(config['T_cross_s']*duration_factor/config['dt_s'])*config['dt_s']
    steps0=math.ceil((T+config['post_pulse_crossings']*config['T_cross_s'])/config['dt_s'])
    every=max(1,math.ceil(steps0/240));steps=math.ceil(steps0/every)*every
    return {**config,'case_id':f'{index:02d}_{pattern}_A{amplitude:.7g}_D{duration_factor:.7g}',
        'pattern_id':pattern,'tau_budget_Pa':amplitude,'duration_s':T,'duration_over_Tcross':T/config['T_cross_s'],
        'nominal_duration_factor':duration_factor,'explicit_steps':steps,'dump_every':every,'map_stage':stage}

def prepare_case(r,g):
    p=OUTPUT/'cases'/r['case_id'];p.mkdir(exist_ok=False);load=prepare_provider(g,r)
    write_typed_data(p/'body.data',g,r);dump(p/'request.json',r);dump(p/'LOAD_PROVIDER.json',load)
    surface=' '.join(map(str,np.flatnonzero(g['surface'])+1));area=float(g['areas'][g['surface']][0]);V=r['dx_m']**3
    text=prelude(g,r,r['dump_every'])+f'''group surface id {surface}
set group surface d_area {area:.17g}
'''+lammps_provider(g,r,load)+f'''fix drag clot viscous {r['rho_kg_m3']*V*r['numerical_drag_rate_s_inv']:.17g}
fix fixed anchors setforce 0 0 0
fix integrate clot nve
timestep {r['dt_s']:.17g}
variable force_mag atom sqrt(v_load_fx^2+v_load_fy^2+v_load_fz^2)
compute force_budget surface reduce sum v_force_mag
variable actual_budget equal c_force_budget/{g['areas'].sum():.17g}
thermo {r['dump_every']*5}
thermo_style custom step time atoms bonds v_amp v_actual_budget ke
thermo_modify lost error lost/bond error flush yes format float %.15g
timer full
'''+dumps(r['dump_every'],r['explicit_steps'],'v_load_fx v_load_fy v_load_fz')
    (p/'in.lammps').write_text(text)
    return dict(path=p,request=r,geometry=g,load=load,steps=r['explicit_steps'])

def analyze_case(case,*,shape_function=spatial_shape):
    start=time.perf_counter();p=case['path'];r=case['request'];g=case['geometry'];n=g['clot_count'];N=len(g['X']);b=g['bonds'];typ=g['bond_types'];M=len(b)
    keys=b[:,0]*N+b[:,1];L0=np.linalg.norm(g['X'][b[:,1]]-g['X'][b[:,0]],axis=1);ev=_events(p/'broken_bonds.local')
    if len(ev):
        ends=np.sort(ev[:,2:4].astype(int)-1,axis=1);ek=ends[:,0]*N+ends[:,1];order=np.argsort(ev[:,1],kind='stable');_,first=np.unique(ek[order],return_index=True);ev=ev[order[first]];ev=ev[np.argsort(ev[:,1],kind='stable')]
        ends=np.sort(ev[:,2:4].astype(int)-1,axis=1);ek=ends[:,0]*N+ends[:,1];idx=np.searchsorted(keys,ek)
        if np.any(idx>=M) or not np.array_equal(keys[idx],ek):raise ValueError('Unknown broken bond')
    else:idx=np.empty(0,int)
    alive=np.ones(M,bool);cursor=0;rows=[];XX=[];VV=[];labels=[];attached_labels=[];step_all=[];peaks={};last_frags=[];provider_error=0.;peak_traction=0.;peak_actual_budget=0.
    with (p/'fragment_components.csv').open('w',newline='') as f:
        fields=['step','time_s','component_id','node_count','volume_m3','mass_kg','attached','com_x_m','com_y_m','com_z_m','vx_m_s','vy_m_s','vz_m_s','speed_m_s'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for step,cols,a in read_dump(p/'trajectory.lammpstrj'):
            a=a[np.argsort(a[:,0])];assert np.array_equal(a[:,0],np.arange(1,N+1)) and np.isfinite(a).all()
            X=a[:,[cols.index(c) for c in ['x','y','z']]];v=a[:,[cols.index(c) for c in ['vx','vy','vz']]]
            assert np.max(abs(X[n:]-g['X'][n:]))<1e-12
            while cursor<len(ev) and ev[cursor,1]<=step:alive[idx[cursor]]=False;cursor+=1
            lab,counts,attached=clot_components(n,b,typ,alive);nc=len(counts)
            ids=np.full(nc,n+1,int);np.minimum.at(ids,lab,np.arange(1,n+1));stable=ids[lab]
            volume=counts*r['dx_m']**3;mass=volume*r['rho_kg_m3'];com=np.stack([np.bincount(lab,weights=X[:n,k])/counts for k in range(3)],axis=1);vc=np.stack([np.bincount(lab,weights=v[:n,k])/counts for k in range(3)],axis=1)
            fr=[]
            for j in range(nc):
                z=dict(step=step,time_s=step*r['dt_s'],component_id=int(ids[j]),node_count=int(counts[j]),volume_m3=float(volume[j]),mass_kg=float(mass[j]),attached=bool(attached[j]),com_x_m=float(com[j,0]),com_y_m=float(com[j,1]),com_z_m=float(com[j,2]),vx_m_s=float(vc[j,0]),vy_m_s=float(vc[j,1]),vz_m_s=float(vc[j,2]),speed_m_s=float(np.linalg.norm(vc[j])));w.writerow(z);fr.append(z)
            metrics=dict(step=step,time_s=step*r['dt_s'],bulk_broken_fraction=float((~alive[typ==1]).mean()),interface_broken_fraction=float((~alive[typ==2]).mean()),interface_survival=float(alive[typ==2].mean()),
                largest_free_fraction=float(counts.max()/n),FI_free=float(1-counts.max()/n),resolved_fragment_count=int((counts>=3).sum()),debris_fraction=float(counts[counts<=2].sum()/n),median_fragment_nodes=float(np.median(counts)),raw_fragment_count=nc,
                attached_clot_volume_fraction=float(counts[attached].sum()/n),max_fragment_speed_m_s=float(np.linalg.norm(vc,axis=1).max()))
            metrics.update(component_metrics(counts,attached,r['dx_m'],r['rho_kg_m3']))
            metrics['failure_mode_detached_v2']=classify_detached_v2(metrics)
            uu=utilization(X,b,typ,alive,L0,r['s00_bulk'],r['s00_interface']);metrics.update(uu);metrics['failure_mode']=classify(metrics)
            for k,value in uu.items():
                if not k.startswith('surviving_') and value is not None:peaks[k]=max(peaks.get(k,0.),value)
            F=a[:,[cols.index('v_load_f'+c) for c in 'xyz']];expected=shape_function(X,r)*g['areas'][:,None]*case['load']['scale_Pa']*float(pulse(step*r['dt_s'],r['duration_s']))
            provider_error=max(provider_error,float(abs(F-expected).max()));budget=float(np.linalg.norm(F,axis=1).sum()/g['areas'].sum());metrics['actual_mean_traction_budget_Pa']=budget;peak_actual_budget=max(peak_actual_budget,budget)
            local=np.linalg.norm(F[g['surface']],axis=1)/g['areas'][g['surface']];peak_traction=max(peak_traction,float(local.max()))
            rows.append(metrics);XX.append(X);VV.append(v);labels.append(stable);attached_labels.append(attached[lab]);step_all.append(step);last_frags=fr
    top=list(read_dump(p/'surviving_bonds.local'))[-1];ends=np.sort(top[2][:,1:3].astype(int)-1,axis=1);actual=np.unique(ends[:,0]*N+ends[:,1]);assert np.array_equal(actual,keys[alive])
    assert alive[typ==3].all(),'Fixed wall topology broke';assert provider_error<max(1e-13,case['load']['initial_integrated_force_budget_N']*1e-10)
    csv_out(p/'fragment_metrics.csv',rows);np.savez_compressed(p/'accepted_trajectory.npz',step=step_all,position_m=XX,velocity_m_s=VV,component_id=labels,attached=attached_labels,clot_count=n,anchors=g['anchors'])
    np.savez_compressed(p/'broken_bonds.npz',events=ev,initial_bond_indices=idx)
    final=rows[-1];result={**final,**peaks,'case_id':r['case_id'],'pattern_id':r['pattern_id'],'map_stage':r['map_stage'],'tau_budget':r['tau_budget_Pa'],
        'peak_local_traction':case['load']['reference_peak_traction_Pa'],'actual_peak_local_traction_Pa':peak_traction,'peak_actual_mean_traction_budget_Pa':peak_actual_budget,
        'duration':r['duration_s'],'duration_over_Tcross':r['duration_over_Tcross'],'nominal_duration_factor':r['nominal_duration_factor'],
        's00_bulk':r['s00_bulk'],'s00_interface':r['s00_interface'],'bulk_target_strength_range_Pa':[4000,8000],'interface_target_strength_range_Pa':[1700,3500],
        'bulk_calibrated_macro_strength_Pa':r['bulk_macro_strength_Pa'],'interface_calibrated_macro_strength_Pa':r['interface_macro_strength_Pa'],
        'clot_nodes':n,'wall_anchors':N-n,'V_free0':n*r['dx_m']**3,'scope':'EXPLORATORY_ONLY','final_fragments':last_frags,
        'utilization_definition':'Primary columns: maximum across saved times of surviving-bond percentile/fraction; final values saved separately',
        'final_utilization':{k:v for k,v in final.items() if k.startswith(('U_','frac_'))},
        'fragment_size_distribution_nodes':[z['node_count'] for z in last_frags],'fragment_volume_distribution_m3':[z['volume_m3'] for z in last_frags],
        'fragment_mass_distribution_kg':[z['mass_kg'] for z in last_frags],'fragment_velocity_statistics':dict(max_speed_m_s=max(z['speed_m_s'] for z in last_frags),median_speed_m_s=float(np.median([z['speed_m_s'] for z in last_frags]))),
        'validation':dict(all_atoms_retained=True,finite_states=True,fixed_wall_anchors=True,wall_bonds_intact=True,native_topology_match=True,provider_max_force_difference_N=provider_error,anchors_excluded=True),
        'runtime':json.loads((p/'PERFORMANCE.json').read_text()),'postprocess_wall_s':time.perf_counter()-start}
    for t,name in [(1,'bulk'),(2,'interface')]:
        sub=ev[typ[idx]==t];result['first_'+name+'_break_time']=float(sub[:,1].min()*r['dt_s']) if len(sub) else None
    dump(p/'result.json',result);return result

def run_one(r,g):
    case=prepare_case(r,g);run_lammps(case);result=analyze_case(case)
    print('MAP_RESULT',r['case_id'],result['failure_mode'],'bulk',result['bulk_broken_fraction'],'interface',result['interface_broken_fraction'],'FI',result['FI_free'],flush=True);return result

def table_row(r):return {k:v for k,v in r.items() if not isinstance(v,(list,dict))}

def run_campaign(config):
    # Explicit gate; cannot quietly use a failed or guessed threshold.
    for kind in ['BULK','INTERFACE']:
        status=json.loads((OUTPUT/(kind+'_CALIBRATION.json')).read_text());assert status['status']=='PASS',kind+' calibration gate'
    assert config['s00_interface']<config['s00_bulk'];g=prepare_interface_geometry(config)
    np.savez_compressed(OUTPUT/'INTERFACE_GEOMETRY.npz',**{k:v for k,v in g.items() if k!='info'});dump(OUTPUT/'INTERFACE_GEOMETRY.json',g['info'])
    start=time.perf_counter();results=[];requests=[];n=0;stop=None
    def batch(reqs):
        with ThreadPoolExecutor(max_workers=config['max_concurrent_cases']) as pool:
            futs={pool.submit(run_one,q,g):q for q in reqs}
            for future in as_completed(futs):
                z=future.result();results.append(z);dump(OUTPUT/'CAMPAIGN_PROGRESS.json',dict(completed=len(results),submitted=len(requests),last_case=z['case_id'],last_mode=z['failure_mode'],elapsed_s=time.perf_counter()-start));dump(OUTPUT/'FAILURE_MAP.json',sorted(results,key=lambda z:z['case_id']));csv_out(OUTPUT/'FAILURE_MAP.csv',[table_row(a) for a in sorted(results,key=lambda z:z['case_id'])])
    amplitudes=list(config['initial_amplitudes_Pa'])
    for amp in amplitudes:
        reqs=[]
        for pat in PATTERNS:
            for duration in config['duration_factors']:
                n+=1;reqs.append(case_request(config,pat,amp,duration,'COARSE',n))
        requests.extend(reqs);csv_out(OUTPUT/'COARSE_MAP_CASES.csv',requests);batch(reqs)
        if amp==10 and all(z['failure_mode']=='PULVERIZATION' for z in results):stop='MODEL_MATERIAL_INCONSISTENCY_AT_10_PA';break
    if stop is None:
        high=[z for z in results if z['tau_budget']==1000]
        if len(high)==9 and all(z['bulk_broken_fraction']==0 and z['interface_broken_fraction']==0 for z in high):
            reqs=[]
            for pat in PATTERNS:
                for duration in config['duration_factors']:
                    n+=1;reqs.append(case_request(config,pat,10000.,duration,'CONDITIONAL_UPPER_BRACKET',n))
            requests.extend(reqs);batch(reqs)
        candidates={p:[] for p in PATTERNS};seen=set()
        for pat in PATTERNS:
            z=[a for a in results if a['pattern_id']==pat]
            for axis,fixed in [('tau_budget','nominal_duration_factor'),('nominal_duration_factor','tau_budget')]:
                for value in sorted(set(a[fixed] for a in z)):
                    line=sorted([a for a in z if a[fixed]==value],key=lambda a:a[axis])
                    for a,b in zip(line,line[1:]):
                        if a['failure_mode']==b['failure_mode']:continue
                        mid=math.sqrt(a[axis]*b[axis]);amp=mid if axis=='tau_budget' else value;dur=value if axis=='tau_budget' else mid;key=(pat,amp,dur)
                        if key in seen:continue
                        seen.add(key);priority=0 if 'NO_FAILURE' in [a['failure_mode'],b['failure_mode']] else 1
                        candidates[pat].append((priority,amp,dur,a['case_id'],b['case_id']))
            candidates[pat].sort()
        selected=[];cap=min(config['max_refinement_cases'],config['max_formal_cases']-n)
        while len(selected)<cap and any(candidates.values()):
            for pat in PATTERNS:
                if candidates[pat] and len(selected)<cap:selected.append((pat,*candidates[pat].pop(0)))
        reqs=[]
        for pat,priority,amp,dur,lower,upper in selected:
            n+=1;q=case_request(config,pat,amp,dur,'BOUNDARY_REFINEMENT',n);q.update(refinement_parent_cases=[lower,upper]);reqs.append(q)
        csv_out(OUTPUT/'ADAPTIVE_REFINEMENT_CASES.csv',reqs);requests.extend(reqs)
        if reqs:batch(reqs)
    else:csv_out(OUTPUT/'ADAPTIVE_REFINEMENT_CASES.csv',[])
    sub=[table_row(z) for z in results if z['failure_mode']=='NO_FAILURE'];csv_out(OUTPUT/'SUBCRITICAL_UTILIZATION.csv',sub)
    dump(OUTPUT/'PERFORMANCE_CAMPAIGN.json',dict(status='COMPLETE' if stop is None else 'STOPPED_INCONSISTENCY',stop_reason=stop,formal_cases=len(results),campaign_wall_s=time.perf_counter()-start,sum_case_wall_s=sum(z['runtime']['wall_seconds'] for z in results),postprocess_sum_wall_s=sum(z['postprocess_wall_s'] for z in results),peak_per_case_RSS_bytes=max(z['runtime']['peak_RSS_bytes'] for z in results),max_concurrent_cases=config['max_concurrent_cases'],MPI_ranks_per_case=config['mpi_ranks'],no_extra_performance_campaign=True))
    return results
