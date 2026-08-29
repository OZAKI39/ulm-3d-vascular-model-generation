-- Isolated adaptive flux validation; no restart and no field export.
simulation_name = 'adaptive_flux_validation'
printRuntimeInfo = true
timing_file = 'musubi_timing.res'
mesh = './mesh/'
scaling = 'diffusive'
logging = { level = 5 }
dx = 1.9999999999999999e-07
dt = 2.4414062499999991e-08
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 23622.32012800001
maximum_iterations = 20
sim_control = {
  time_control = { max = {iter=maximum_iterations}, interval={iter=1} },
  abort_criteria = { stop_file = 'stop' }
}
physics = { dt=dt, rho0=rho0_phy }
identify = { label='adaptive_flux_test', kind='fluid', layout='d3q19', relaxation='bgk' }
fluid = { kinematic_viscosity=nu_phy, bulk_viscosity=bulk_viscosity_phy }
initial_condition = { pressure=pressure_reference_phy, velocityX=0.0, velocityY=0.0, velocityZ=0.0 }
boundary_condition = {
  {label='inlet', kind='adaptive_flux_pressure', mass_flowrate=8.124344950169123e-13},
  {label='outlet', kind='pressure_eq', pressure=pressure_reference_phy},
  {label='wall_x0', kind='wall'}, {label='wall_x1', kind='wall'},
  {label='wall_y0', kind='wall'}, {label='wall_y1', kind='wall'}
}
