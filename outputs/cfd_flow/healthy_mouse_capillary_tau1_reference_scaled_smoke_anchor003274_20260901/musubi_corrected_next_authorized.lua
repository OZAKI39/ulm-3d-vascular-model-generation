-- Fresh Tau1 reference-pressure smoke; numerical offset is not physiological pressure.
simulation_name = 'tau1_reference_scaled_smoke'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh/'
scaling = 'diffusive'
logging = {level=5}
maximum_iterations = 5000
dx = 1.9999999999999999e-07
dt = 2.0387359836901121e-09
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 3387510.7199999993
function outlet_01_pressure(x,y,z,t) return 3387525.2649781005 end
function outlet_02_pressure(x,y,z,t) return 3387642.9245492225 end
function outlet_03_pressure(x,y,z,t) return 3387497.0193733261 end
physics = {dt=dt, rho0=rho0_phy}
identify = {label='ROI003274_tau1', kind='fluid', layout='d3q19', relaxation='bgk'}
fluid = {kinematic_viscosity=nu_phy, bulk_viscosity=bulk_viscosity_phy}
initial_condition = {pressure=pressure_reference_phy, velocityX=0.0, velocityY=0.0, velocityZ=0.0}
boundary_condition = {
  {label='wall', kind='wall_libb'},
  {label='inlet', kind='adaptive_flux_pressure', mass_flowrate=2.8901803804796421e-12},
  {label='outlet_01', kind='pressure_eq', pressure=outlet_01_pressure},
  {label='outlet_02', kind='pressure_eq', pressure=outlet_02_pressure},
  {label='outlet_03', kind='pressure_eq', pressure=outlet_03_pressure}
}
sim_control = {
  time_control={max={iter=maximum_iterations}, interval={iter=100}},
  abort_criteria={stop_file='stop'}
}
tracking = {
  {label='state_000000', folder='tracking/state_000000/', variable={'density_phy','velocity_phy'}, shape={kind='all'}, time_control={min={iter=0}, max={iter=0}, interval={iter=1}}, output={format='asciiSpatial', use_get_point=false}},
  {label='state_000001', folder='tracking/state_000001/', variable={'density_phy','velocity_phy'}, shape={kind='all'}, time_control={min={iter=1}, max={iter=1}, interval={iter=1}}, output={format='asciiSpatial', use_get_point=false}},
  {label='state_000010', folder='tracking/state_000010/', variable={'density_phy','velocity_phy'}, shape={kind='all'}, time_control={min={iter=10}, max={iter=10}, interval={iter=1}}, output={format='asciiSpatial', use_get_point=false}},
  {label='state_000100', folder='tracking/state_000100/', variable={'density_phy','velocity_phy'}, shape={kind='all'}, time_control={min={iter=100}, max={iter=100}, interval={iter=1}}, output={format='asciiSpatial', use_get_point=false}},
  {label='state_001000', folder='tracking/state_001000/', variable={'density_phy','velocity_phy'}, shape={kind='all'}, time_control={min={iter=1000}, max={iter=1000}, interval={iter=1}}, output={format='asciiSpatial', use_get_point=false}},
  {label='state_005000', folder='tracking/state_005000/', variable={'density_phy','velocity_phy'}, shape={kind='all'}, time_control={min={iter=5000}, max={iter=5000}, interval={iter=1}}, output={format='asciiSpatial', use_get_point=false}}
}
restart = {write='restart/', timeformat={use_iter=true},
  time_control={min={iter=4999}, max={iter=maximum_iterations}, interval={iter=1}}
}
