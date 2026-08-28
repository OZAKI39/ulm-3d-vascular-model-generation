-- Generated production configuration; official Musubi syntax.
simulation_name = 'roi003274_steady_lbm'
printRuntimeInfo = true
timing_file = 'mus_timing.res'
mesh = '../../musubi_recovery_anchor003274_20260828_162530/seeder/mesh/'
scaling = 'diffusive'
logging = { level = 5 }

dx = 1.9999999999999999e-07
dt = 2.4414062499999991e-08
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = (2.0 / 3.0) * nu_phy
pressure_reference_phy = 23622.32012800001
maximum_iterations = 1000000

-- Requested upstream parabola, retained explicitly for reproducibility.
-- The effective official BC is mfr_eq so exact Q takes priority.
inlet_center = { 0.0001036737050000898, 5.7554233238768548e-05, 0.00015880639894086851 }
inlet_inward_normal = { -0.025172709087112898, -0.67208421514314431, -0.74004671641230813 }
inlet_equivalent_radius = 1.5776895780401971e-06
inlet_maximum_velocity = 0.00019677116015072346
function requested_parabolic_velocity(x, y, z, t)
  local rx = x - inlet_center[1]
  local ry = y - inlet_center[2]
  local rz = z - inlet_center[3]
  local axial = rx*inlet_inward_normal[1] + ry*inlet_inward_normal[2] + rz*inlet_inward_normal[3]
  rx = rx - axial*inlet_inward_normal[1]
  ry = ry - axial*inlet_inward_normal[2]
  rz = rz - axial*inlet_inward_normal[3]
  local factor = math.max(0.0, 1.0 - (rx*rx + ry*ry + rz*rz)/(inlet_equivalent_radius^2))
  return { inlet_inward_normal[1]*inlet_maximum_velocity*factor,
            inlet_inward_normal[2]*inlet_maximum_velocity*factor,
            inlet_inward_normal[3]*inlet_maximum_velocity*factor }
end

function outlet_01_pressure(x, y, z, t) return 23636.865106101286 end
function outlet_02_pressure(x, y, z, t) return 23754.524677223188 end
function outlet_03_pressure(x, y, z, t) return 23608.619501326699 end

sim_control = {
  time_control = {
    max = { iter = maximum_iterations, clock = 3600 },
    interval = { iter = 100 }
  },
  abort_criteria = {
    stop_file = 'stop',
    steady_state = true,
    convergence = {
      variable = { 'pressure_phy', 'vel_mag_phy' },
      shape = { kind = 'all' },
      reduction = { 'average', 'average' },
      time_control = { min = { iter = 0 }, max = { iter = maximum_iterations }, interval = { iter = 100 } },
      norm = 'average',
      nvals = 100,
      absolute = true,
      condition = {
        { threshold = 0.001, operator = '<=' },
        { threshold = 1.0000000000000001e-09, operator = '<=' }
      }
    }
  }
}

physics = { dt = dt, rho0 = rho0_phy }
identify = { label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }
fluid = {
  kinematic_viscosity = nu_phy,
  bulk_viscosity = bulk_viscosity_phy
}
initial_condition = { pressure = pressure_reference_phy, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }

boundary_condition = {
  { label = 'wall', kind = 'wall_libb' },
  { label = 'inlet', kind = 'mfr_eq', mass_flowrate = 8.124344950169123e-13 },

  { label = 'outlet_01', kind = 'pressure_eq', pressure = outlet_01_pressure },
  { label = 'outlet_02', kind = 'pressure_eq', pressure = outlet_02_pressure },
  { label = 'outlet_03', kind = 'pressure_eq', pressure = outlet_03_pressure }
}

restart = { write = 'restart/' }
