"""CPU regression tests. All manufactured samples/processes are SYNTHETIC,
never included in scientific calibration results or GPU budget usage.
"""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT,config_read,output,write_json,read_json,fingerprint
from py_scripts.fluid_physics.units import Units,KB,kelvin,thermal_density_units,consistency,pressure_target_star,viscosity_only_alternative
from py_scripts.fluid_physics.legacy_case import freeze_case,temporal_residual,time_mean,flow_closure
from py_scripts.fluid_physics.analysis import fit_viscosity,periodic_shape,block_statistics,correlation_time_samples,native_pressure,read_csv
from py_scripts.fluid_physics.calibration import package_units,verify_package
from py_scripts.fluid_physics.runner import bounded_process,group_members,run_attempt,ledger_read,used_budget
from py_scripts.fluid_physics.boundary_plan import build_boundary_plan


class UnitsTests(unittest.TestCase):
    def setUp(self):self.u=thermal_density_units(5e-7,1056,8,1,1,298.15)
    def test_all_quantities_roundtrip(self):
        for name in self.u.scales:
            with self.subTest(name=name):np.testing.assert_allclose(self.u.to_star(self.u.to_si([.1,1,7],name),name),[.1,1,7],rtol=1e-13)
    def test_density(self):self.assertAlmostEqual(self.u.to_si(8,'mass_density'),1056)
    def test_viscosities_and_flows(self):
        rho,nu,Q=1056,3.27e-6,2.7369132390905703e-15
        rs=self.u.to_star(rho,'mass_density');ns=self.u.to_star(nu,'kinematic_viscosity');qs=self.u.to_star(Q,'volume_flow')
        self.assertAlmostEqual(self.u.to_si(rs*ns,'dynamic_viscosity'),rho*nu)
        self.assertAlmostEqual(self.u.to_si(rs*qs,'mass_flow')/(rho*Q),1)
    def test_celsius_kelvin_energy(self):
        self.assertEqual(kelvin(25),298.15);self.assertAlmostEqual(self.u.scales['energy']/(KB*298.15),1)
    def test_nonpositive_or_nonfinite(self):
        for bad in (0,-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):Units(bad,1,1)
    def test_reject_incompatible_thermal_viscous(self):
        d=consistency(self.u,n_star=8,m_star=1,kBT_star=1,rho_si=1056,temperature_K=298.15,nu_star=1.5,nu_si=3.27e-6)
        self.assertTrue(d['thermal_consistent']);self.assertFalse(d['viscosity_consistent']);self.assertFalse(d['all_consistent'])
        alt=viscosity_only_alternative(self.u,1.5,3.27e-6,1)
        self.assertGreater(alt['implied_temperature_K'],1e6);self.assertFalse(alt['selection_allowed'])
    def test_nonunit_particle_mass(self):
        u=thermal_density_units(5e-7,1056,4,2,1,298.15)
        self.assertAlmostEqual(u.to_si(4*2,'mass_density'),1056)
    def test_negative_gauge_and_reference(self):
        p=pressure_target_star(150,-13.700626673311461,0,self.u)
        self.assertLess(p,150);self.assertAlmostEqual(self.u.to_si(p-150,'pressure'),-13.700626673311461)
        self.assertAlmostEqual(pressure_target_star(150,10,-10,self.u)-150,self.u.to_star(20,'pressure'))


class LegacyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=config_read(PROJECT_ROOT/'py_scripts/fluid_physics.yaml');cls.case,cls.a,cls.h,cls.g=freeze_case(cls.c)
    def test_real_source_chain(self):
        self.assertEqual(self.case['accepted_iteration'],598755);self.assertTrue(self.case['source_matching']['generated_lua_equals_accepted'])
        self.assertNotEqual(self.case['source_matching']['old_binary_sha256'],self.case['source_matching']['current_config_binary_sha256'])
        self.assertTrue(self.case['source_matching']['accepted_flux_diagnostics_recomputed'])
    def test_actual_legacy_denominators(self):
        fraction=next(r for r in self.a['criteria'] if r['name']=='flow_fraction_drift')
        density=next(r for r in self.a['criteria'] if r['name']=='Q_density_consistency')
        self.assertIn('sum Qout',fraction['formula']);self.assertIn('/abs(Q_rho_u_over_rho0)',density['formula'])
    def test_original_cut_cap_measurement_distinct(self):
        for p in self.case['ports']:
            original=np.array(p['original_physical_cut_reference']['center_m']);cap=np.array(p['patch']['center_m']);measured=np.array(p['measurement_planes']['central']['origin_m'])
            self.assertGreater(np.linalg.norm(original-cap),1e-5);self.assertLess(np.linalg.norm(cap-measured),2e-6)
    def test_source_conflict_fails(self):
        c=copy.deepcopy(self.c);c['geometry_manifest_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'SOURCE_CONFLICT'):freeze_case(c)
    def test_mass_volume_source(self):
        self.assertAlmostEqual(self.case['target_mass_flow_kg_s']/(self.case['density_kg_m3']*self.case['target_volume_flow_m3_s']),1)
    def test_port_identity_and_planes(self):
        self.assertEqual(len({p['patch']['port_id'] for p in self.case['ports']}),4)
        for p in self.case['ports']:
            self.assertGreater(np.dot(p['patch']['outward_normal'],p['measurement_planes']['central']['unit_normal']),.99999)
            self.assertGreater(np.linalg.norm(np.array(p['patch']['center_m'])-p['measurement_planes']['central']['origin_m']),1e-6)
    def test_rigid_shape_and_unchanged(self):
        self.assertLess(self.case['source_matching']['rotated_wall_triangle_max_error_m'],5e-11);self.assertFalse(self.case['source_matching']['geometry_modified'])
    def test_source_numeric_pressure(self):
        self.assertGreater(self.case['pressure_reference']['legacy_numerical_offset_pa'],3e6)
        self.assertIsNone(self.case['pressure_reference']['thermodynamic_absolute_pressure_pa'])
        self.assertLess(self.case['outlet_gauge_pressures_pa']['outlet_03'],0)
    def test_temporal_not_target_error(self):
        self.assertEqual(temporal_residual(100,100,1),0);self.assertGreater(abs(100-10)/10,1)
        r=next(x for x in self.a['criteria'] if x['name']=='R_pressure');self.assertIn('NOT_TARGET_ERROR',r['category'])
    def test_insufficient_physical_window(self):
        with self.assertRaisesRegex(ValueError,'WINDOW_INSUFFICIENT'):time_mean([0,.1],[1,1],.2)
        self.assertAlmostEqual(time_mean([0,1,3],[0,2,2],3),5/3)
    def test_signed_backflow_preserved(self):
        self.assertEqual(flow_closure(1,[.9,.12,-.02])['significant_backflow'],[False,False,False])
        self.assertTrue(flow_closure(1,[1.2,-.2,0])['significant_backflow'][1])
    def test_boundary_finite_and_explicit_gaps(self):
        u,_,_=package_units(self.case,self.c);b=build_boundary_plan(self.case,self.c,u)
        self.assertFalse(b['cpp_cuda_modified']);self.assertEqual(b['status'],'REQUIRES_EXTENSION')
        self.assertEqual(len(b['ports']),4)
        for p in b['ports']:
            self.assertEqual(p['state'],'REQUIRES_EXTENSION');self.assertGreater(p['measurement_to_cap_axial_distance_m'],1e-6)
    def test_locked_thermal_scale_and_bubble_ratios(self):
        u,m,cs=package_units(self.case,self.c)
        self.assertGreater(m['required_nu_star'],400);np.testing.assert_allclose(m['bubble_to_spacing_ratio'],[8,16]);self.assertEqual(len(cs['candidates']),2)


class StatisticsTests(unittest.TestCase):
    def test_analytic_recovery_nonunit_mass(self):
        y=np.arange(.125,8,.25);nu=2.7;F=.2;m=3.4
        v=F/(m*nu)*periodic_shape(y,8,.25)+.123
        r=fit_viscosity(y,v,F,m,8,.25);self.assertAlmostEqual(r['nu_star'],nu,places=12);self.assertLess(r['relative_rms'],1e-12)
        self.assertLess(v[4]-.123,0);self.assertGreater(v[-5]-.123,0)
    def test_numerator_force_not_acceleration(self):
        y=np.arange(.125,8,.25);v=.1/(2*4)*periodic_shape(y,8,.25)
        self.assertAlmostEqual(fit_viscosity(y,v,.1,2,8,.25)['nu_star'],4)
        self.assertAlmostEqual(fit_viscosity(y,v,.1,1,8,.25)['nu_star'],8)
    def test_correlated_noise_effective_samples(self):
        rng=np.random.default_rng(123);x=np.zeros(20000)
        for i in range(1,len(x)):x[i]=.99*x[i-1]+rng.normal()
        r=block_statistics(np.arange(len(x),dtype=float),x,min_duration=1)
        self.assertGreater(r['correlation_time_samples'],80);self.assertLess(r['effective_sample_estimate'],len(x)/80)
        self.assertGreater(r['block_samples'],400)
    def test_insufficient_and_unstable_series(self):
        r=block_statistics(np.arange(100,dtype=float),np.arange(100,dtype=float),min_duration=20)
        self.assertEqual(r['status'],'WINDOW_INSUFFICIENT')
        self.assertEqual(block_statistics([0,1],[1,1],min_duration=1)['status'],'WINDOW_INSUFFICIENT')
    def test_nonfinite_cannot_fit(self):
        with self.assertRaises(ValueError):fit_viscosity([0,1,2,3],[1,2,float('nan'),4],1,1,8)
    def test_pressure_kinetic_volume_nonunit_mass(self):
        s=np.array([(0,2,2,0,0,100)],dtype=[('time',float),('kBT',float),('vx',float),('vy',float),('vz',float),('num_particles',float)])
        v=np.array([(0,30)],dtype=[('time',float),('pressure',float)])
        p,t=native_pressure(s,v,mass=2,volume=10)
        self.assertAlmostEqual(p[0],3+10*(2-2/3));self.assertAlmostEqual(t[0],(2-2/3)*100/99)
        s['time']=.002
        p2,_=native_pressure(s,v,2,10,dt=.002);self.assertEqual(p[0],p2[0])
        v['time']=1
        with self.assertRaisesRegex(ValueError,'TIMESTAMP'):native_pressure(s,v,2,10)


class ProtectionBudgetTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='synthetic_cpu_',dir=PROJECT_ROOT/'test_code/outputs/fluid_physics');self.d=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def test_source_symlink_rejected(self):
        p=self.d/'link';p.symlink_to(PROJECT_ROOT.parent/'ulm_3D_vascular',target_is_directory=True)
        with self.assertRaises(Exception):output(p/'must_not_write.json')
    def test_exclusive_json_and_finite(self):
        p=self.d/'a.json';write_json(p,{'a':1});self.assertEqual(read_json(p),{'a':1})
        with self.assertRaises(FileExistsError):write_json(p,{'a':2})
        with self.assertRaises(ValueError):write_json(self.d/'nan.json',{'a':float('nan')})
    def test_missing_file_fails(self):
        with self.assertRaisesRegex(ValueError,'INPUT_MISSING'):read_csv(self.d/'missing.csv')
    def test_export_manifest_roundtrip_and_tamper(self):
        from py_scripts.fluid_physics.common import sha256_file
        p=self.d/'payload.json';write_json(p,{'category':'SYNTHETIC_CPU_ONLY','value':2})
        write_json(self.d/'package_sha256.json',{'payload.json':sha256_file(p)})
        verify_package(self.d);p.write_text('{"value":3}')
        with self.assertRaisesRegex(ValueError,'HASH_MISMATCH'):verify_package(self.d)
    def test_timeout_only_own_group(self):
        unrelated=subprocess.Popen([sys.executable,'-c','import time; time.sleep(20)'],start_new_session=True)
        d=self.d/'timeout';d.mkdir()
        code="import subprocess,sys,time,pathlib; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(20)']); pathlib.Path('child_pid').write_text(str(p.pid)); time.sleep(20)"
        try:
            r=bounded_process([sys.executable,'-c',code],d,.6,{'graceful_margin_s':.2},monitor_gpu=False)
            self.assertTrue(r['timeout']);self.assertFalse(r['remaining_own_group_processes']);self.assertIsNone(unrelated.poll())
            self.assertLess(r['elapsed_monotonic_s'],1.2)
            self.assertLess(r['started_at_utc'],r['ended_at_utc'])
        finally:os.killpg(unrelated.pid,15);unrelated.wait()
    def test_persistent_cache_retry_budget(self):
        c={'campaign_id':'synthetic','budget':{'task_limit_s':2,'campaign_limit_s':2,'minimum_launch_budget_s':.1,'graceful_margin_s':.1}}
        campaign=self.d/'campaign'
        make=lambda d:[sys.executable,'-c','pass']
        d,r,hit=run_attempt(campaign,c,'a',{'m':1,'units':1,'binary':'a','precision':'single','script':'a'},make,.8,monitor_gpu=False)
        self.assertFalse(hit);used=used_budget(ledger_read(campaign,'synthetic',c['budget']))
        d2,r2,hit=run_attempt(campaign,c,'a',{'m':1,'units':1,'binary':'a','precision':'single','script':'a'},make,.8,monitor_gpu=False)
        self.assertTrue(hit);self.assertEqual(used,used_budget(ledger_read(campaign,'synthetic',c['budget'])))
        for k,v in [('m',2),('units',2),('binary','b'),('precision','double'),('script','b')]:
            ident={'m':1,'units':1,'binary':'a','precision':'single','script':'a'};ident[k]=v
            with self.assertRaisesRegex(RuntimeError,'EXISTING_DIFFERENT'):run_attempt(campaign,c,'a',ident,make,.8,monitor_gpu=False)
        run_attempt(campaign,c,'a',{'m':2},make,.8,retry_failed=True,monitor_gpu=False)
        self.assertGreater(used_budget(ledger_read(campaign,'synthetic',c['budget'])),used)
    def test_crashed_reservation_consumes_budget(self):
        ledger={'attempts':[{'reserved_s':600,'status':'RUNNING'},{'reserved_s':600,'charged_s':3,'status':'FAILED'}]}
        self.assertEqual(used_budget(ledger),603)
    def test_budget_above_authorization_rejected(self):
        c=config_read(PROJECT_ROOT/'py_scripts/fluid_physics.yaml');c['budget']['campaign_limit_s']=3601
        import yaml
        p=self.d/'bad.yaml';p.write_text(yaml.safe_dump(c))
        with self.assertRaisesRegex(ValueError,'budget'):config_read(p)
    def test_help_does_not_import_mirheo(self):
        for module in ('py_scripts.prepare_fluid_physics','py_scripts.calibrate_dpd_fluid','test_code.review_fluid_physics'):
            p=subprocess.run([sys.executable,'-B','-m',module,'--help'],cwd=PROJECT_ROOT,capture_output=True,text=True,timeout=10)
            self.assertEqual(p.returncode,0);self.assertNotIn('Mirheo version',p.stdout)


class EosReportTests(unittest.TestCase):
    """Manufactured report inputs test errors/coverage, never scientific outputs."""
    def test_minimal_eos_coverage_and_json(self):
        from py_scripts.fluid_physics.reporting import eos_analysis
        ts=[]
        for name,n in [('eos_low',7.6),('equilibrium',8.0),('eos_high',8.4)]:
            ts.append({'candidate_id':'thermal_baseline','kind':'equilibrium','task_id':name,'actual_n_star':n,'parameters':{'m_star':1},'pressure':{'mean':10*n,'ci95_halfwidth':.1,'status':'SUFFICIENT'},'temperature_status':'PASS_PROPOSED'})
        u=Units(1,1,1);case={'outlet_gauge_pressures_pa':{'a':1,'b':100,'c':-200}}
        r=eos_analysis(ts,case,u,{})
        self.assertAlmostEqual(r['fit']['dp_drho_star'],10);self.assertEqual(r['target_pressure_coverage_status'],'OUTSIDE_MEASURED_EOS_RANGE')
        self.assertIsNone(r['density_target_solution']);json.dumps(r,allow_nan=False)
        self.assertEqual(eos_analysis(ts[:2],case,u,{})['status'],'INCONCLUSIVE')
    def test_sparse_eos_does_not_gain_precision(self):
        from py_scripts.fluid_physics.reporting import eos_analysis
        ts=[{'candidate_id':'thermal_baseline','kind':'equilibrium','task_id':name,'actual_n_star':n,'parameters':{'m_star':1},'pressure':{'mean':10*n,'ci95_halfwidth':30,'status':'WINDOW_INSUFFICIENT'},'temperature_status':'INCONCLUSIVE'} for name,n in [('eos_low',7.6),('equilibrium',8),('eos_high',8.4)]]
        r=eos_analysis(ts,{'outlet_gauge_pressures_pa':{'a':1}},Units(1,1,1),{})
        self.assertEqual(r['status'],'INCONCLUSIVE')
    def test_sensitivity_near_mean_but_broad_ci_fails(self):
        from py_scripts.fluid_physics.reporting import sensitivity
        tasks=[{'task_id':n,'viscosity':{'nu_star':1,'sampling_status':'SUFFICIENT','ci95_star':[.5,1.5]}} for n in ('flow','flow_half_dt','flow_half_force')]
        r=sensitivity(tasks,{'sensitivity_relative_difference':.1});self.assertNotEqual(r['status'],'PASS_PROPOSED')


if __name__=='__main__':unittest.main()
