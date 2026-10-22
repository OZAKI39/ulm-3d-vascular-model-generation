"""Official APECSS refinement; no changes to old tables or radial equations."""
import time, dataclasses, concurrent.futures, numpy as np, h5py
from mevo_mb.bubbles.apecss import ApecssBubbleResponseModel, APECSS_COMMIT
from scripts.straight_tube_exp1.f10ad_multifrequency.build import case
from scripts.straight_tube_exp1.f10s.common import context
from .common import *
def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-80))
def solve(arg):
    p,level,tag=arg;cfg=read(CONFIG)['low_pressure']
    c=dataclasses.replace(case(2.5e6,p),cycles=cfg['refinement_cycles'][level],samples_per_cycle=cfg['samples_per_cycle'][level],tolerance=cfg['tolerances'][level])
    folder=BASE/'low_pressure'/tag
    if (folder/'RESULT.json').exists():
        md=read(folder/'RESULT.json')
        h=np.load(folder/'cycles.npz')['harmonics_cosine_m3'] if md['status']=='COMPLETE' else None
    else:
        try:md,h,cyc=ApecssBubbleResponseModel(ROOT/'third_party/build/ultrasound_apecss').solve(c,folder)
        except Exception:
            md=read(folder/'RESULT.json');h=None
    return dict(pressure_Pa=p,level=level,status=md['status'],error=md.get('error'),harmonics=h,metrics=md.get('metrics'),case=md['case'],wall_s=md['elapsed_s'],path=str(folder.relative_to(BASE)))
def main():
    verify();start=time.perf_counter();cfg=read(CONFIG)['low_pressure'];old=context()['response']
    audit=dict(historical_response=old.metadata,actual_fixed_interval_min_Pa=float(old.knots[0]),actual_fixed_interval_max_Pa=float(old.knots[-1]),new_frequency_tables=[])
    for f in [2400000,2500000,2600000]:
        for shell in [.8,1.,1.2]:
            p=AD/str(f)/f'response_shell_{shell:.1f}.h5'
            with h5py.File(p) as h: audit['new_frequency_tables'].append(dict(path=str(p),sha256=sha(p),pressure_range_Pa=[float(h['pressure_Pa'][0]),float(h['pressure_Pa'][-1])],metadata=__import__('json').loads(h.attrs['metadata_json'])))
    failures=[]
    for r in read(AC/'ALL_RUNS.json'):
        if r['status']=='PASS':continue
        entry=dict(id=r['id'],reason=r.get('reasons'),pressure_scope='Independent CPU reconstruction at LAST ACCEPTED state; NOT original failed kernel pressure')
        try:
            with h5py.File(AC/r['trajectory_path']) as h:x=h['positions_m'][-1]
            c=context();p,g=c['basis'].query_basis(x);u=c['basis'].vector(r['task']['amplitudes'],r['task']['phases']);pressure=u@p
            entry['reconstructed_min_Pa']=float(abs(pressure).min())
            try:old.query(pressure);entry['fail_closed']=False
            except ValueError as e:entry.update(fail_closed=True,error=str(e))
        except Exception as e:entry['audit_error']=repr(e)
        failures.append(entry)
    audit['historical_failures']=failures;write(BASE/'LOW_PRESSURE_TABLE_AUDIT.json',audit)
    nodes=cfg['nodes_Pa'];args=[(p,k,f'node_{i:02}_level{k}') for i,p in enumerate(nodes) for k in range(3)]
    mids=[(a+b)/2 for a,b in zip(nodes[:-1],nodes[1:])]
    args += [(p,2,f'mid_{i:02}') for i,p in enumerate(mids)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(solve,args))
    write(BASE/'LOW_PRESSURE_RAW.json',rows)
    checks=[];omega=2*np.pi*2.5e6;harmonics=np.arange(1,9)
    for i,p in enumerate(nodes):
        rr=[r for r in rows if r['pressure_Pa']==p];h=[np.asarray(r['harmonics']) for r in rr]
        if any(r['status']!='COMPLETE' for r in rr):
            checks.append(dict(pressure_Pa=p,radial_checks=False,reason='APECSS refinement failed; no extrapolation'));continue
        # Force-relevant comparisons: complex fundamental F1 and all derivative harmonics F2.
        energy=[float(np.sum(abs(omega*harmonics*z)**2)) for z in h]
        checks.append(dict(pressure_Pa=p,coarse_fine_complex=relative(h[0],h[1]),fine_finest_complex=relative(h[1],h[2]),fundamental_relative=relative(h[1][:1],h[2][:1]),F2_energy_relative=abs(energy[1]-energy[2])/max(energy[2],1e-80),phase_rad=float(np.angle(h[2][0])),volume_amplitude_m3=float(abs(h[2][0])),volume_susceptibility=h[2][0]/p,derivative_energy_m6_s2=energy[2],radial_checks=all(r['metrics']['valid'] and r['metrics']['harmonic_energy_resolved'] for r in rr)))
    interp=[]
    for i,p in enumerate(mids):
        needed=[r for r in rows if r['level']==2 and r['pressure_Pa'] in [p,nodes[i],nodes[i+1]]]
        if any(r['status']!='COMPLETE' for r in needed):
            interp.append(dict(pressure_Pa=p,lower_Pa=nodes[i],upper_Pa=nodes[i+1],valid=False,reason='Missing refined endpoint'));continue
        ref=np.asarray(next(r['harmonics'] for r in rows if r['pressure_Pa']==p));left=np.asarray(next(r['harmonics'] for r in rows if r['pressure_Pa']==nodes[i] and r['level']==2));right=np.asarray(next(r['harmonics'] for r in rows if r['pressure_Pa']==nodes[i+1] and r['level']==2));pred=(left+right)/2
        interp.append(dict(pressure_Pa=p,lower_Pa=nodes[i],upper_Pa=nodes[i+1],valid=True,complex_relative=relative(pred,ref),F2_energy_relative=abs(np.sum(abs(harmonics*pred)**2)-np.sum(abs(harmonics*ref)**2))/np.sum(abs(harmonics*ref)**2)))
    junction=np.asarray(next(r['harmonics'] for r in rows if r['pressure_Pa']==.1 and r['level']==2));oldh=old.query(np.array([.1+0j]))[0]
    junction_error=relative(junction,oldh)
    # Admit only a contiguous independently converged suffix; failed smaller nodes stay illegal.
    def node_ok(r):return r['radial_checks'] and r.get('fine_finest_complex',1)>0 and r['fine_finest_complex']<cfg['complex_response_relative_limit'] and r['fundamental_relative']<cfg['force_relative_limit'] and r['F2_energy_relative']<cfg['force_relative_limit']
    legal_start=len(nodes)
    for i in range(len(nodes)-1,-1,-1):
        interval_ok=i==len(nodes)-1 or (interp[i]['valid'] and interp[i]['complex_relative']<cfg['interpolation_relative_limit'] and interp[i]['F2_energy_relative']<cfg['derivative_energy_relative_limit'])
        if not node_ok(checks[i]) or not interval_ok:break
        legal_start=i
    valid_nodes=nodes[legal_start:]
    converged=force=ip=len(valid_nodes)>=2 and valid_nodes[0]<.1
    gates=dict(LOW_PRESSURE_APECSS_CONVERGED=converged,LOW_PRESSURE_FORCE_VALIDATED=force,LOW_PRESSURE_INTERPOLATION_VALIDATED=ip,HISTORICAL_TABLE_UNCHANGED=sha(OLD/'inputs/response.h5')==read(BASE/'PREREGISTRATION.json')['protected_sha256'][str((OLD/'inputs/response.h5').relative_to(ROOT))],HISTORICAL_JUNCTION_PASS=junction_error<cfg['historical_junction_relative_limit'],HISTORICAL_2P5MHZ_REGRESSION_PASS='PENDING_BASELINE',LOW_PRESSURE_DYNAMIC_TEST_PASS='PENDING_DIAGNOSTIC')
    # Standalone diagnostic table, never substituted into production controls.
    meta=dict(frequency_Hz=2500000.,radius_m=R,validated=converged and force and ip,scope='LOW_PRESSURE_DIAGNOSTIC_ONLY; not an admitted production extension',identity=identity())
    with h5py.File(BASE/'low_pressure/diagnostic_response.h5','x') as h:
        h['pressure_Pa']=valid_nodes;h['harmonics_m3']=np.array([next(r['harmonics'] for r in rows if r['pressure_Pa']==p and r['level']==2) for p in valid_nodes]);h.attrs['metadata_json']=__import__('json').dumps(meta)
    write(BASE/'LOW_PRESSURE_ANALYSIS.json',dict(diagnostic_pressure_range_Pa=[valid_nodes[0],valid_nodes[-1]] if valid_nodes else None,raw_solver_sum_s=sum(r['wall_s'] for r in rows),failed_radial_solves=sum(r['status']!='COMPLETE' for r in rows),checks=checks,interpolation=interp,junction_relative_error=junction_error,gates=gates,production_extension_admitted=False,decision='Retain original fail-closed domain pending complete dynamic/regression gates; no source kernels or historical table edited',APECSS_commit=APECSS_COMMIT,wall_s=time.perf_counter()-start,identity=identity()))
    print('LOW_PRESSURE',gates,'junction',junction_error,flush=True)
if __name__=='__main__':main()
