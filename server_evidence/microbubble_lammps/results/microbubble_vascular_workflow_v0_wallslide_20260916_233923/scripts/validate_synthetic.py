from pathlib import Path
import csv,json,numpy as np
S=Path(__file__).resolve().parents[1]
assert 'PLANE_A_B_C_D_E_PASS' in (S/'validation/SYNTHETIC_NATIVE.log').read_text()
rows=list(csv.DictReader((S/'validation/SYNTHETIC_CURVED.csv').open()));results=[]
for shape in ['cylinder','sphere']:
    series=[]
    for dt in [1e-4,5e-5,2.5e-5]:
        r=[x for x in rows if x['shape']==shape and float(x['dt_max'])==dt];assert r and abs(float(r[-1]['time'])-.01)<1e-14
        assert min(float(x['gap']) for x in r)>=1e-10-2e-14
        assert all(0<float(x['dt'])<=dt for x in r)
        last=r[-1];x=np.array([float(last[k]) for k in 'xyz']);rc=5e-6-1e-6-1e-10;angle=1e-4*.01/rc;ref=rc*np.array([np.cos(angle),np.sin(angle),0.]);err=float(np.linalg.norm(x-ref))
        series.append(dict(dt_max=dt,steps=len(r),trajectory_error_m=err,max_correction_m=float(last['max_correction']),total_correction_m=float(last['total_correction']),constraint_events=int(last['constraint_events']),path_length_m=float(last['path_length']),retry_count=int(last['retries']),min_gap_m=min(float(x['gap']) for x in r)))
    for key in ['trajectory_error_m','max_correction_m','total_correction_m']:
        assert series[2][key]<series[1][key]<series[0][key],(shape,key,series)
    results.append(dict(shape=shape,status='PASS',refinement=series))
out=dict(status='PASS',flat_A_to_E='PASS',plane_tangent_and_outward_relative_tolerance=1e-12,curved=results,position_correction_qualification='maximum per-step and cumulative corrections decrease under dt/2 and dt/4',native_test_source='src/test_wall_constraint.cpp')
(S/'validation/SYNTHETIC_WALL_VALIDATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
