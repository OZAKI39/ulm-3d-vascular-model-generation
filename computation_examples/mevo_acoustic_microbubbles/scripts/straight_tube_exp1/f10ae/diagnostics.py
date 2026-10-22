"""Supplemental read-only pressure/gradient/force and low-pressure scaling audit."""
import time,numpy as np,h5py
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import context
from .common import *
def main():
    begin=time.perf_counter();rows=[]
    for r in read(BASE/'ROBUSTNESS_RUNS.json'):
        if r['status']!='PASS':continue
        t=r['task'];c=context(t.get('frequency',2.5e6),t.get('shell',1.));u=c['basis'].vector(t['amplitudes'],t['phases'])
        with h5py.File(BASE/r['trajectory_path']) as h:
            x=h['positions_m'][0];f1=h['force_samples/primary_force_N'][0];f2=h['force_samples/secondary_force_N'][0];actual=h['force_samples/pressure_Pa'][0]
        pp,gg=c['basis'].query_basis(x);grad=np.einsum('i,inj->nj',u,gg);p=u@pp;mean=f1.mean(0);norm=float(np.linalg.norm(mean))
        rows.append(dict(id=r['id'],frequency=t.get('frequency',2.5e6),shell=t.get('shell',1.),offset_m=t.get('offset_m',[0,0,0]),pressure_min_Pa=float(abs(actual).min()),pressure_max_Pa=float(abs(actual).max()),pressure_cpu_recorded_relative=float(np.linalg.norm(p-actual)/max(np.linalg.norm(actual),1e-80)),gradient_mean_complex_Pa_m=grad.mean(0),gradient_mean_norm_Pa_m=float(np.linalg.norm(grad,axis=1).mean()),primary_mean_N=mean,primary_direction=mean/norm if norm>0 else None,primary_mean_magnitude_N=float(np.linalg.norm(f1,axis=1).mean()),secondary_mean_magnitude_N=float(np.linalg.norm(f2,axis=1).mean()),shape=r['targets'],scope='Pressure and forces from saved first step; gradient reconstructed from the exact matching cached field'))
    write(BASE/'ROBUSTNESS_FIELD_FORCE.json',dict(rows=rows,wall_s=time.perf_counter()-begin,identity=identity()))
    low=read(BASE/'LOW_PRESSURE_ANALYSIS.json');v=[r for r in low['checks'] if r.get('radial_checks') and r['pressure_Pa']>=.001 and r['pressure_Pa']<=.1];p=np.array([r['pressure_Pa'] for r in v]);vol=np.array([r['volume_amplitude_m3'] for r in v]);energy=np.array([r['derivative_energy_m6_s2'] for r in v])
    write(BASE/'ZERO_LIMIT_DIAGNOSTIC.json',dict(pressure_range_Pa=[float(p.min()),float(p.max())],volume_amplitude_loglog_slope=float(np.polyfit(np.log(p),np.log(vol),1)[0]),derivative_energy_loglog_slope=float(np.polyfit(np.log(p),np.log(energy),1)[0]),volume_secant_relative_spread=float(np.ptp(vol/p)/np.mean(vol/p)),scope='Empirical approach toward zero over validated nonzero pressures; no extrapolation below .001 Pa. OFF exact zero remains analytically distinct from nonzero controls near a node.',identity=identity()))
if __name__=='__main__':main()
