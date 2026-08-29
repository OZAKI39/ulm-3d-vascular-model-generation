-- Isolated adaptive flux validation; no restart and no field export.
simulation_name = 'adaptive_flux_validation'
printRuntimeInfo = true
timing_file = 'musubi_timing.res'
mesh = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_ideal_plane_inlet_preflight_anchor003274_20260829_120444/seeder/mesh/'
scaling = 'diffusive'
logging = { level = 5 }
dx = 1.9999999999999999e-07
dt = 2.4414062499999991e-08
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 23622.32012800001
maximum_iterations = 30
sim_control = {
  time_control = { max = {iter=maximum_iterations}, interval={iter=1} },
  abort_criteria = { stop_file = 'stop' }
}
physics = { dt=dt, rho0=rho0_phy }
identify = { label='adaptive_flux_test', kind='fluid', layout='d3q19', relaxation='bgk' }
fluid = { kinematic_viscosity=nu_phy, bulk_viscosity=bulk_viscosity_phy }
initial_condition = { pressure=pressure_reference_phy, velocityX=0.0, velocityY=0.0, velocityZ=0.0 }
function outlet_01_pressure(x,y,z,t) return 23636.865106101286 end
function outlet_02_pressure(x,y,z,t) return 23754.524677223188 end
function outlet_03_pressure(x,y,z,t) return 23608.619501326699 end
boundary_condition = {
  {label='wall', kind='wall_libb'},
  {label='inlet', kind='adaptive_flux_pressure', mass_flowrate=8.124344950169123e-13},
  {label='outlet_01', kind='pressure_eq', pressure=outlet_01_pressure},
  {label='outlet_02', kind='pressure_eq', pressure=outlet_02_pressure},
  {label='outlet_03', kind='pressure_eq', pressure=outlet_03_pressure}
}
