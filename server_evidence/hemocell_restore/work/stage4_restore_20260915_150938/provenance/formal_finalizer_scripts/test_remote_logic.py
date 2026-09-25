"""Synthetic mathematics/IO tests only. NO lattice initialization or timesteps."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from remote_common import PORTS,OFFSETS,FIELD_DTYPE,field_path
from convergence import exact_window,evaluate_point,initial_state,field_residual
from verify_run import reconstruct_flux

ROOT=Path(__file__).resolve().parent
CONTRACT_PATH=(ROOT.parent/'frozen_contracts/convergence_contract.json') if (ROOT.parent/'frozen_contracts').exists() else Path('/home/lzy/projects/compre_output/step3b/20260913_023318/contracts/convergence_contract.json')
C=json.loads(CONTRACT_PATH.read_text())

def observations(end=700000):
    names=['iteration','time_s','rho_min','rho_max','rho_mean','u_max_m_s','Mach_max','mean_speed_m_s','total_mass','control_volume_mass','relative_mass_drift','BOUNDARY_PROFILE_Q_TARGET']
    names += [f'{p}_{k}dx'+('_dx2' if r else '') for p in PORTS for k in OFFSETS for r in [0,1]]
    names += [f'mass_outward_g{g}' for g in range(24)]
    names += ['instantaneous_flow_closure','control_volume_mass_balance','ghost_rho_min','ghost_rho_max','ghost_Mach_max']
    data=np.zeros(end//100+1,dtype=[(n,'f8') for n in names]);data['iteration']=np.arange(len(data))*100;data['time_s']=data['iteration']*C['dt_s']
    for name in ['rho_min','rho_max','rho_mean','ghost_rho_min','ghost_rho_max']:data[name]=1
    data['mean_speed_m_s']=1e-4;data['total_mass']=2e-12;data['control_volume_mass']=1.5e-12
    data['BOUNDARY_PROFILE_Q_TARGET']=C['Qtarget_m3_s'];data['control_volume_mass_balance'][0]=np.nan
    fractions=[1,.1,.7,.2]
    for p,pname in enumerate(PORTS):
        for k,offset in enumerate(OFFSETS):
            for r in [0,1]:
                data[f'{pname}_{offset}dx'+('_dx2' if r else '')]=C['Qtarget_m3_s']*fractions[p]
                data[f'mass_outward_g{6*p+2*k+r}']=(-1 if p==0 else 1)*C['rho_kg_m3']*C['Qtarget_m3_s']*fractions[p]
    return data

class FrozenLogicTests(unittest.TestCase):
    def test_exact_linear_window_non_grid_endpoints(self):
        t=np.arange(0,5,.2);y=3*t+2
        self.assertAlmostEqual(exact_window(t,y,.31,3.77),3*(.31+3.77)/2+2,places=13)
        with self.assertRaises(ValueError):exact_window(t,y,-.01,2)

    def test_early_evaluation_and_full_confirmation_interval(self):
        d=observations();s=initial_state()
        row,s,details=evaluate_point(d[d['iteration']<=235000],235000,C,s,None)
        self.assertFalse(row['eligible']);self.assertFalse(details['stop'])
        for step in range(240000,360001,5000):
            row,s,details=evaluate_point(d[d['iteration']<=step],step,C,s,(0.,0.,0.))
            self.assertTrue(row['all_gates_pass'])
            self.assertEqual(details['stop'],step==360000)
        self.assertEqual(s['candidate_iteration'],240000);self.assertEqual(s['confirmation_iteration'],360000)

    def test_any_intermediate_failure_resets_confirmation(self):
        d=observations();s=initial_state()
        for step in range(240000,425001,5000):
            fv=(.2,0.,0.) if step==300000 else (0.,0.,0.)
            row,s,detail=evaluate_point(d[d['iteration']<=step],step,C,s,fv)
            if step==300000:self.assertIsNone(s['candidate_iteration']);self.assertEqual(s['consecutive_pass_count'],0)
            self.assertEqual(detail['stop'],step==425000)
        self.assertEqual(s['first_converged_iteration'],240000);self.assertEqual(s['candidate_iteration'],305000)

    def test_backflow_is_not_hidden_by_absolute_values(self):
        d=observations(240000)
        for k in OFFSETS:
            d[f'Qout01_{k}dx']=-.001*C['Qtarget_m3_s']
            d[f'Qout01_{k}dx_dx2']=d[f'Qout01_{k}dx']
        row,state,detail=evaluate_point(d,240000,C,initial_state(),(0.,0.,0.))
        self.assertFalse(detail['gates']['outlet_01_mean_direction']);self.assertFalse(row['all_gates_pass'])

    def test_plane_and_quadrature_fail_independently(self):
        d=observations(240000);d['Qin_2dx']*=1.04;d['Qin_2dx_dx2']*=1.04
        row,_,detail=evaluate_point(d,240000,C,initial_state(),(0.,0.,0.))
        self.assertFalse(detail['gates']['cross_plane_flux_spread']);self.assertTrue(detail['gates']['quadrature_error'])
        d=observations(240000);d['Qin_4dx_dx2']*=.98
        _,_,detail=evaluate_point(d,240000,C,initial_state(),(0.,0.,0.))
        self.assertTrue(detail['gates']['cross_plane_flux_spread']);self.assertFalse(detail['gates']['quadrature_error'])

    def test_transient_storage_can_balance_without_steady_closure(self):
        d=observations(240000)
        for g in range(6,24):d[f'mass_outward_g{g}']*=.9
        for p in PORTS[1:]:
            for k in OFFSETS:
                d[f'{p}_{k}dx']*=.9;d[f'{p}_{k}dx_dx2']*=.9
        d['control_volume_mass']+=.1*C['rho_kg_m3']*C['Qtarget_m3_s']*d['time_s']
        row,_,detail=evaluate_point(d,240000,C,initial_state(),(0.,0.,0.))
        self.assertLess(row['R_CV'],1e-9);self.assertGreater(row['R_flow'],.09)

    def test_missing_sample_or_nonfinite_data_rejected(self):
        d=observations(240000)
        with self.assertRaises(ValueError):evaluate_point(np.delete(d,37),240000,C,initial_state(),(0.,0.,0.))
        d['Qin_4dx'][123]=np.nan
        with self.assertRaises(ValueError):evaluate_point(d,240000,C,initial_state(),(0.,0.,0.))

    def test_full_field_change_detected_despite_identical_mean_speed(self):
        with tempfile.TemporaryDirectory(prefix='step3b_synthetic_fields_') as td:
            run=Path(td);(run/'contracts').mkdir();n=182694
            a=np.zeros(n,dtype=FIELD_DTYPE);a['index']=np.arange(n);a['rho']=1;a['u'][:,0]=1e-6
            b=a.copy();b['u'][:n//2,0]=0;b['u'][:n//2,1]=1e-6
            np.savez(run/'contracts/control_volume_nodes.npz',flat_index=a['index'],control_volume_mask=np.ones(n,'u1'))
            for step,data in [(0,a),(C['formulas']['field_delta_steps'],b)]:
                p=field_path(run,step);p.parent.mkdir(parents=True,exist_ok=True)
                with p.open('wb') as f:np.array([n],dtype='<u8').tofile(f);data.tofile(f)
            ru,rp,source_ru=field_residual(run,C['formulas']['field_delta_steps'],C)
            self.assertGreater(ru,.5);self.assertEqual(source_ru,0.);self.assertEqual(rp,0.)

    def test_independent_trilinear_mass_product(self):
        points=np.array([[x,y,z] for z in [0,1] for y in [0,1] for x in [0,1]])
        sampled=np.zeros(8,dtype=FIELD_DTYPE);sampled['index']=np.arange(8);sampled['rho']=1+.1*points[:,2];sampled['u']=points
        q=np.array([[2,.25,.5,.75,2.,0.,0.,1.]])
        flux,mass=reconstruct_flux(q,sampled,np.array([2,2,2]),.2,.1,1056.)
        expected=2*.75*.2**3/.1
        self.assertAlmostEqual(flux[2],expected,places=15)
        self.assertAlmostEqual(mass[2],1056*1.075*expected,places=12)

    def test_horizon_does_not_force_pass(self):
        d=observations();d['Qin_4dx']*=.8
        row,_,detail=evaluate_point(d,700000,C,initial_state(),(0.,0.,0.))
        self.assertFalse(detail['stop']);self.assertFalse(row['all_gates_pass']);self.assertFalse(detail['gates']['R_inlet'])

if __name__=='__main__':unittest.main(verbosity=2)
