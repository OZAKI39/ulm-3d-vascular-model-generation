-- Frozen-geometry, force-driven periodic Poiseuille validation.
simulation_name = 'periodic_pipe_force_axis_n16'
printRuntimeInfo = true
timing_file = 'tracking/musubi_timing.res'
mesh = 'mesh/'
scaling = 'diffusive'
logging = { level = 5 }
maximum_iterations = 2000
dx = 1.9721119725502463e-07
dt = 2.373794941574722e-08
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 23622.32012800001
radius = 1.5776895780401971e-06
target_mean = 0.00035
axis = { 0, 0, 1 }
center = { 2.5243033248643153e-05, 2.5243033248643153e-05, 2.5243033248643153e-05 }

function radial_squared(x,y,z)
  local rx=x-center[1]; local ry=y-center[2]; local rz=z-center[3]
  local axial=rx*axis[1]+ry*axis[2]+rz*axis[3]
  rx=rx-axial*axis[1]; ry=ry-axial*axis[2]; rz=rz-axial*axis[3]
  return rx*rx+ry*ry+rz*rz
end
function vel_analy(x,y,z,t)
  return 2.0*target_mean*math.max(0.0,1.0-radial_squared(x,y,z)/(radius*radius))
end
function initial_velocity_x(x,y,z,t) return axis[1]*vel_analy(x,y,z,t) end
function initial_velocity_y(x,y,z,t) return axis[2]*vel_analy(x,y,z,t) end
function initial_velocity_z(x,y,z,t) return axis[3]*vel_analy(x,y,z,t) end

physics = { dt = dt, rho0 = rho0_phy }
identify = { label = 'periodic_pipe_force', kind = 'fluid_incompressible', layout = 'd3q19', relaxation = 'bgk' }
fluid = { kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }
initial_condition = {
  pressure = pressure_reference_phy,
  velocityX = initial_velocity_x, velocityY = initial_velocity_y, velocityZ = initial_velocity_z
}
boundary_condition = { { label = 'wall', kind = 'wall_libb' } }
glob_source = { force = { 0, 0, 3884423.6432636492 }, force_order = 2 }

variable = {
  { name = 'vel_analy', ncomponents = 1, vartype = 'st_fun', st_fun = vel_analy },
  { name = 'vel_error', ncomponents = 1, vartype = 'operation',
     operation = { kind = 'difference', input_varname = {'vel_mag_phy','vel_analy'} } }
}

sim_control = {
  time_control = { max = { iter = maximum_iterations }, interval = { iter = 50 } },
  abort_criteria = { stop_file = 'stop' }
}

tracking = {
  { label = 'mean_velocity', folder = 'tracking/', variable = {'velocity_phy'},
     shape = { kind = 'all' }, reduction = {'average'},
     time_control = { min = {iter=0}, max = {iter=maximum_iterations}, interval = {iter=50} },
     output = { format = 'ascii' } },
  { label = 'profile', folder = 'tracking/', variable = {'vel_analy','vel_error'},
     shape = { kind = 'all' }, reduction = {'l2norm','l2norm'},
     time_control = { min = {iter=0}, max = {iter=maximum_iterations}, interval = {iter=50} },
     output = { format = 'ascii' } },
  { label = 'safety', folder = 'tracking/', variable = {'vel_mag','pdf'},
     shape = { kind = 'all' }, reduction = {'max','min'},
     time_control = { min = {iter=0}, max = {iter=maximum_iterations}, interval = {iter=50} },
     output = { format = 'ascii' } }
}
