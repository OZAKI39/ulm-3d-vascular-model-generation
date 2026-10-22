"""Independent returned-file, numerical-domain and force-derivative review."""
import time
import subprocess
import unittest
import numpy as np
import h5py
from .common import *
from .diagnostics import field,primary,jacobian
from .test_admission import Coverage

def main():
 start=time.perf_counter();checks={};evidence={};cfg=read(CONFIG);c=context();p=policy()
 test=unittest.TestResult();unittest.defaultTestLoader.loadTestsFromTestCase(Coverage).run(test)
 checks['duration_regression_and_incomplete_window_rejection']=test.wasSuccessful()
 evidence['test_failures']=test.failures+test.errors
 checks['unchanged_physics_src_and_prior_configs']=not subprocess.check_output(['git','diff',cfg['baseline_commit'],'--','src','configs/straight_tube_exp1/f10q','configs/straight_tube_exp1/f10r','scripts/straight_tube_exp1/f10q','scripts/straight_tube_exp1/f10r'],cwd=ROOT,text=True)
 manifests=read(BASE/'REMOTE_MANIFEST.json')['files']
 checks['all_returned_files_verified']=all(sha(BASE/k)==v['sha256'] for k,v in manifests.items())
 for cp in BASE.glob('*/CAMPAIGN.json'):
  for rp in (cp.parent/'runs').glob('*/RESULT.json'):
   r=read(rp);h=h5py.File(rp.parent/'trajectory.h5');x=h['positions_m'][:];t=h['time_s'][:];rows=h['step_metrics'][:]
   i,j=np.triu_indices(32,1);gap=np.linalg.norm(x[:,i]-x[:,j],axis=-1)-2*p['radius_m'];coord=(x-c['basis'].origin)/c['basis'].spacing
   passed=bool(r['status']=='PASS' and len(x)-1==r['accepted_steps'] and rows[:,0].all() and abs(t[-1]-r['requested_s'])<1e-13 and np.min(gap)>=-p['hard_geometry_m'] and np.array_equal(x[0],c['x0']) and not rows[:,5].any() and np.all(coord>=1) and np.all(coord<np.array(c['basis'].shape)-2))
   for col,limit in [(1,p['mobility_residual_m']),(2,p['contact_overlap_m']),(3,p['contact_correction_m']),(4,p['contact_correction_m']),(12,p['primal_violation_m']),(13,p['dual_negative_N']),(15,p['complementarity_equivalent_gap_m']),(17,p['actual_geometry_complementarity_m'])]:passed=passed and rows[:,col].max()<=limit
   checks[r['name']+'_original_basic_thresholds']=bool(passed);h.close()
 # Diagnostic quadrature and derivatives checked independently of stored cubic caches.
 u=np.array([0,0,1036.2016581041792,1036.2016581041792*np.exp(1j*np.deg2rad(6))]);pts=np.array([[0,0,0],[5e-6,1e-6,-1e-6]])
 a,b=c['discs'].evaluate(pts,48);aa,bb=c['discs'].evaluate(pts,96)
 evidence['quadrature_max_P_basis_difference']=float(abs(a-aa).max());evidence['quadrature_max_gradient_basis_difference']=float(abs(b-bb).max())
 checks['disc_quadrature']=bool(np.max(abs(a-aa))<1e-12 and np.max(abs(b-bb))<1e-7)
 point=np.array([3e-6,1e-6,-2e-6]);h=1e-9;P,G=field(c,[point],u)
 fd=np.array([(field(c,[point+np.eye(3)[i]*h],u)[0][0]-field(c,[point-np.eye(3)[i]*h],u)[0][0])/(2*h) for i in range(3)])
 relative=float(np.linalg.norm(fd-G[0])/np.linalg.norm(G[0]));checks['analytic_pressure_gradient_consistency']=relative<1e-7;evidence['gradient_FD_relative_error']=relative
 data=read(BASE/'diagnostics/TRAPS.json');checks['trap_J_halved_difference']=all(r['force_J_halfstep_max_difference_N_m']<abs(np.array(r['force_J_eigenvalues_N_m'])).min()*.001 for r in data['records'])
 phase=[r for r in read(BASE/'diagnostics/Y_SCOPE.json')['rows'] if r['kind']=='phase']
 checks['single_source_global_phase_invariance']=all(np.allclose(r['primary_force_N'],phase[0]['primary_force_N'],rtol=1e-10,atol=1e-24) for r in phase)
 checks['no_off_frequency_APECSS_forces']=all(r['primary_force_N'] is None for r in read(BASE/'diagnostics/Y_SCOPE.json')['rows'] if r['frequency_Hz']!=2.5e6)
 off=R/'pilot_v1/runs/OFF_DT1/trajectory.h5'
 with h5py.File(off) as f: checks['historical_OFF_exactly_stationary']=bool(np.all(f['positions_m'][:]==f['positions_m'][0]))
 write(BASE/'analysis/INDEPENDENT_REVIEW.json',dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,evidence=evidence,wall_s=time.perf_counter()-start,identity=identity()))
 print(checks);assert all(checks.values()),checks
if __name__=='__main__':main()
