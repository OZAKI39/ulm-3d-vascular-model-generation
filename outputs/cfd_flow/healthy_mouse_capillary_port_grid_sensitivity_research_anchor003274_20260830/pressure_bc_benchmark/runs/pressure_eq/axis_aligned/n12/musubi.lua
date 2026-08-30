-- Isolated fixed-geometry pressure-BC benchmark; not production CFD.
simulation_name = 'pipe_pressure_eq_12'
printRuntimeInfo = true
timing_file = 'timing.res'
mesh = 'mesh/'
scaling = 'diffusive'
logging = { level = 5 }
maximum_iterations = 8000
dx = 2.6294826300669951e-07
dt = 4.2200798961328387e-08
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 23622.32012800001
function outlet_pressure(x,y,z,t) return 23622.32012800001 end
physics = { dt = dt, rho0 = rho0_phy }
identify = { label = 'pressure_bc_pipe', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }
fluid = { kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }
initial_condition = { pressure = 23659.090616192032, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }
boundary_condition = {
  { label = 'wall', kind = 'wall_libb' },
  { label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = 2.8901803804796421e-12 },
  { label = 'outlet', kind = 'pressure_eq', pressure = outlet_pressure }
}
sim_control = {
  time_control = { max = { iter = maximum_iterations }, interval = { iter = 100 } },
  abort_criteria = { stop_file = 'stop' }
}
