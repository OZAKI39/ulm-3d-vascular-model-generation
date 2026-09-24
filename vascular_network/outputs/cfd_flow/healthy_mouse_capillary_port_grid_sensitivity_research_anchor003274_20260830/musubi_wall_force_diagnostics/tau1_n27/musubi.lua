-- N27 isolation: identical physical Pipe Force problem, tau exactly 1.
simulation_name = 'periodic_pipe_force_axis_n27_tau1'
printRuntimeInfo = true
timing_file = 'tracking/musubi_timing.res'
mesh = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_port_grid_sensitivity_research_anchor003274_20260830/dimensionless_kernel/periodic_pipe_force/cases/axis_n27/mesh/'
scaling = 'diffusive'
logging = {level=5}
maximum_iterations = 5000
dx = 1.1686589466964423e-07
dt = 6.9610791727504479e-10
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = 2.1799999999999999e-06
pressure_reference_phy = 23622.32012800001
radius = 1.5776895780401971e-06
target_mean = 0.00035
center = {1.4958834517714462e-05,1.4958834517714462e-05,1.4958834517714462e-05}
function radial_squared(x,y,z)
  local rx=x-center[1]; local ry=y-center[2]
  return rx*rx+ry*ry
end
function vel_analy(x,y,z,t)
  return 2.0*target_mean*math.max(0.0,1.0-radial_squared(x,y,z)/(radius*radius))
end
function initial_velocity_z(x,y,z,t) return vel_analy(x,y,z,t) end
physics = {dt=dt, rho0=rho0_phy}
identify = {label='periodic_pipe_force_tau1', kind='fluid_incompressible', layout='d3q19', relaxation='bgk'}
fluid = {kinematic_viscosity=nu_phy, bulk_viscosity=bulk_viscosity_phy}
initial_condition = {pressure=pressure_reference_phy, velocityX=0.0, velocityY=0.0, velocityZ=initial_velocity_z}
boundary_condition = {{label='wall', kind='wall_libb'}}
glob_source = {force={0, 0, 3884423.6432636492}, force_order=2}
variable = {
  {name='vel_analy', ncomponents=1, vartype='st_fun', st_fun=vel_analy},
  {name='vel_error', ncomponents=1, vartype='operation',
    operation={kind='difference', input_varname={'vel_mag_phy','vel_analy'}}}
}
sim_control = {
  time_control={max={iter=maximum_iterations}, interval={iter=50}},
  abort_criteria={
    stop_file='stop', steady_state=true,
    convergence={
      variable={'vel_mag_phy','vel_error'}, shape={kind='all'},
      reduction={'average','l2norm'},
      time_control={min={iter=200}, max={iter=maximum_iterations}, interval={iter=50}},
      norm='average', nvals=5, absolute=false,
      condition={{threshold=1.0e-4,operator='<='},{threshold=1.0e-4,operator='<='}}
    }
  }
}
tracking = {
  {label='mean_velocity', folder='tracking/', variable={'velocity_phy'},
    shape={kind='all'}, reduction={'average'},
    time_control={min={iter=0}, max={iter=maximum_iterations}, interval={iter=50}},
    output={format='ascii'}},
  {label='profile', folder='tracking/', variable={'vel_analy','vel_error'},
    shape={kind='all'}, reduction={'l2norm','l2norm'},
    time_control={min={iter=0}, max={iter=maximum_iterations}, interval={iter=50}},
    output={format='ascii'}},
  {label='safety', folder='tracking/', variable={'vel_mag','pdf'},
    shape={kind='all'}, reduction={'max','min'},
    time_control={min={iter=0}, max={iter=maximum_iterations}, interval={iter=50}},
    output={format='ascii'}},
  {label='cross_section', folder='tracking/', variable={'velocity_phy'},
    shape={kind='canoND', object={
      origin={1.3147413150334976e-05,1.3147413150334976e-05,1.4958834517714462e-05},
      vec={{3.6228427347589708e-06,0.0,0.0},{0.0,3.6228427347589708e-06,0.0}}
    }},
    time_control={min={iter=0}, max={iter=maximum_iterations}, interval={iter=50}},
    output={format='asciiSpatial', use_get_point=false}}
}
