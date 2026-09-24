-- Fresh repaired Fine logical simulation; first hard stop is 5000.
simulation_name = 'tau1_reference_scaled_fine_postfix'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/fine/seeder/mesh/'
scaling = 'diffusive'
logging = {level=5}
maximum_iterations = 5000
dx = 1.5384615384615385e-07
dt = 1.2063526530710723e-09
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 5724893.1168
function outlet_01_pressure(x,y,z,t) return 5724907.6617781017 end
function outlet_02_pressure(x,y,z,t) return 5725025.3213492231 end
function outlet_03_pressure(x,y,z,t) return 5724879.4161733268 end
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
  time_control={max={iter=maximum_iterations}, interval={iter=2023}},
  abort_criteria={stop_file='/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix/stop'}
}
restart = {write='/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix/restart/',
  timeformat={use_iter=true},
  time_control={min={iter=5000}, max={iter=maximum_iterations}, interval={iter=5000}}
}
