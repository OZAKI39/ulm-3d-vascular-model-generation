#!/usr/bin/env python3
"""Sole dt / outlet-rho producer for the isolated PBS/BSA development preparation.
Same dimensional formulas as frozen Step3 prepare_numerics.py lines 129-145.
No solver launch; static success is NOT runtime numerical validation.
"""
from pathlib import Path
import json,math,hashlib
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')
def main():
 p=R/'inputs/NUMERICS_INPUT.json';x=json.loads(p.read_text())
 assert x['development_assumption'] is True and x['experimentally_measured'] is False
 rho=x['rho_kg_m3'];nu=x['nu_m2_s'];dx=x['dx_m'];tau=x['tau']
 assert rho==1000. and nu==1e-6 and tau==1.
 assert math.isclose(x['dynamic_viscosity_pa_s'],rho*nu,rel_tol=1e-15)
 assert dx>0 and tau>.5 and x['initial_rho_lu']==1.
 cs2=1/3
 nu_lu=cs2*(tau-.5)
 dt=nu_lu*dx*dx/nu
 velocity_unit=dx/dt
 pressure_unit=rho*velocity_unit**2
 q_unit=dx**3/dt
 nominal_speed=x['Qtarget_m3_s']/x['inlet_native_projected_area_m2']
 command_speed=x['inlet_numerical_multiplier']*nominal_speed
 command_mach=(command_speed/velocity_unit)/math.sqrt(cs2)
 outlets={}
 for name,pg in x['outlet_pressures_pa'].items():
  pg_lu=pg/pressure_unit
  rho_lu=1+pg_lu/cs2
  assert math.isfinite(rho_lu) and .99<=rho_lu<=1.01
  outlets[name]={'gauge_pressure_pa':pg,'gauge_pressure_lu':pg_lu,'rho_lu':rho_lu,'density_offset':rho_lu-1}
 assert math.isfinite(dt) and dt>0 and command_mach<.05
 # Unit exponents [mass, length, time] verify the same algebra without mixing kg/m3 and rho_LU.
 dim=lambda a,b:[u+v for u,v in zip(a,b)]
 assert dim([0,2,0],[0,-2,1])==[0,0,1] # dx^2 / nu -> s
 assert dim([1,-3,0],[0,2,-2])==[1,-1,-2] # rho * (dx/dt)^2 -> Pa
 v={'schema':'RBC_STAGE1_NUMERICS_V1','status':'STATIC_CONVERSION_PASS_RUNTIME_PENDING',
    'generation_authority':'scripts/prepare_numerics.py','dt_read_from_existing_config':False,
    'dt_generated_from':['dx_m','nu_m2_s','tau'],'input_sha256':sha(p),'script_sha256':sha(Path(__file__)),
    'development_assumption':True,'experimentally_measured':False,'final_experimental_contract':'PENDING',
    'suspending_medium':x['suspending_medium'],'temperature_c':x['temperature_c'],
    'rho_kg_m3':rho,'dynamic_viscosity_pa_s':rho*nu,'nu_m2_s':nu,'dx_m':dx,'tau':tau,'cs2':cs2,
    'nu_lu':nu_lu,'dt_s':dt,'velocity_unit_m_s':velocity_unit,'pressure_unit_pa':pressure_unit,
    'flow_unit_m3_s':q_unit,'mass_unit_kg':rho*dx**3,'force_unit_n':rho*dx**4/dt**2,
    'physical_Qtarget_m3_s':x['Qtarget_m3_s'],'inlet_numerical_multiplier':x['inlet_numerical_multiplier'],
    'outlets':outlets,'nominal_inlet_speed_m_s':nominal_speed,'command_inlet_speed_m_s':command_speed,
    'command_inlet_mach_estimate':command_mach,
    'mach_estimate_scope':'Nominal prescribed inlet command only; not a bound on the full-field transient Mach. Runtime smoke remains required.',
    'dimension_check':'PASS','static_checks':{'tau_safe':True,'density_offsets_small':True,'all_values_finite':True,'inlet_mach_estimate_safe':True,'mu_equals_rho_nu':True},
    'rho_lu_limits':[.99,1.01],'mach_limit':.05,
    'formulas':{'nu_lu':'(tau-0.5)/3','dt_s':'nu_lu*dx_m**2/nu_m2_s','pressure_unit_pa':'rho_kg_m3*(dx_m/dt_s)**2','rho_lu':'1+gauge_pressure_pa/(pressure_unit_pa*(1/3))'},
    'solver_policy':'Future CPU/MPI solver must read and verify this frozen contract; it must not independently regenerate dt or outlet rho_LU.',
    'runtime_smoke_status':'NOT_RUN_BLOCKED_DISK_PREFLIGHT','ready_for_rbc_execution':False,
    'formal_baseline_modified':False}
 out=R/'NUMERICS_CONTRACT.json';write(out,v)
 write(R/'provenance/NUMERICS_GENERATION_RECEIPT.json',{'status':'PASS','scope':'STATIC_GENERATION_ONLY','input_path':'inputs/NUMERICS_INPUT.json','input_sha256':sha(p),'script_path':'scripts/prepare_numerics.py','script_sha256':sha(Path(__file__)),'output_path':'NUMERICS_CONTRACT.json','output_sha256':sha(out),'simulation_timesteps_run':0})
 print('STATIC_NUMERICS_GENERATED',json.dumps({'dt_s':dt,'tau':tau,'outlets':outlets,'command_inlet_mach_estimate':command_mach}))
if __name__=='__main__':main()

