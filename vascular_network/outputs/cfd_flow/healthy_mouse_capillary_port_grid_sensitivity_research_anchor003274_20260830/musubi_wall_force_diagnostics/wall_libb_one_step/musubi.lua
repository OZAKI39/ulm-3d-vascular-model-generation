-- One-step arbitrary-q wall_libb compiled-binary oracle.
simulation_name = 'musubi_wall_libb_one_step_oracle'
printRuntimeInfo = true
timing_file = 'musubi_timing.res'
mesh = 'mesh/'
scaling = 'diffusive'
logging = {level=5}
maximum_iterations = 1
L = 4.6746357867857693e-07
velocity_scale = 14.019487932141525
function initial_velocity_x(x,y,z,t) return velocity_scale*(0.0011 + 0.0003*x/L + 0.0002*y/L - 0.0001*z/L) end
function initial_velocity_y(x,y,z,t) return velocity_scale*(-0.0007 + 0.0001*x/L - 0.00025*y/L + 0.00015*z/L) end
function initial_velocity_z(x,y,z,t) return velocity_scale*(0.0004 - 0.0002*x/L + 0.00015*y/L + 0.0003*z/L) end
physics = {dt=8.3359602886574595e-09, rho0=1056}
identify = {label='wall_oracle', kind='fluid', layout='d3q19', relaxation='bgk'}
fluid = {kinematic_viscosity=3.27e-06, bulk_viscosity=2.1799999999999999e-06}
initial_condition = {
  pressure=23622.32012800001,
  velocityX=initial_velocity_x, velocityY=initial_velocity_y, velocityZ=initial_velocity_z
}
boundary_condition = {{label='wall', kind='wall_libb'}}
sim_control = {
  time_control={max={iter=maximum_iterations}, interval={iter=1}},
  abort_criteria={stop_file='stop'}
}
restart = {
  write='restart/',
  time_control={min={iter=0}, max={iter=1}, interval={iter=1}}
}
