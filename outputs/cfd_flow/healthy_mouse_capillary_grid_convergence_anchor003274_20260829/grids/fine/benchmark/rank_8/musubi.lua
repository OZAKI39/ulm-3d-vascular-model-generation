-- Isolated healthy adaptive-flux grid convergence run.
simulation_name = 'gc_f_b8'
printRuntimeInfo = true
timing_file = 'timing.res'
mesh = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_grid_convergence_anchor003274_20260829/grids/fine/seeder/mesh/'
scaling = 'diffusive'
logging = { level = 5 }
maximum_iterations = 2000
dx = 1.5384615384615385e-07
dt = 1.4446190828402362e-08
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 23622.32012800001

function outlet_01_pressure(x,y,z,t) return 23636.865106101286 end
function outlet_02_pressure(x,y,z,t) return 23754.524677223188 end
function outlet_03_pressure(x,y,z,t) return 23608.619501326699 end

physics = { dt = dt, rho0 = rho0_phy }
identify = { label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }
fluid = { kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }
initial_condition = { pressure = pressure_reference_phy, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }
boundary_condition = {
  { label = 'wall', kind = 'wall_libb' },
  { label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = 2.8901803804796421e-12 },
  { label = 'outlet_01', kind = 'pressure_eq', pressure = outlet_01_pressure },
  { label = 'outlet_02', kind = 'pressure_eq', pressure = outlet_02_pressure },
  { label = 'outlet_03', kind = 'pressure_eq', pressure = outlet_03_pressure }
}
sim_control = {
  time_control = { max = { iter = maximum_iterations }, interval = { iter = 100 } },
  abort_criteria = { stop_file = 'stop' }
}

