-- One-step compiled force-only oracle at production density/scaling.
simulation_name = 'musubi_force_one_step_oracle'
printRuntimeInfo = true
timing_file = 'tracking/musubi_timing.res'
mesh = 'mesh/'
scaling = 'diffusive'
logging = { level = 5 }
maximum_iterations = 1
physics = { dt = 8.3359602886574595e-09, rho0 = 1056 }
identify = {label='force_oracle', kind='fluid_incompressible', layout='d3q19', relaxation='bgk'}
fluid = {kinematic_viscosity=3.27e-06, bulk_viscosity=2.1799999999999999e-06}
initial_condition = {pressure=23622.32012800001, velocityX=0.0, velocityY=0.0, velocityZ=0.0}
glob_source = {force={0, 0, 3884423.6432636492}, force_order=2}
sim_control = {
  time_control={max={iter=maximum_iterations}, interval={iter=1}},
  abort_criteria={stop_file='stop'}
}
tracking = {
  {label='momentum', folder='tracking/', variable={'velocity_phy','density_phy'},
    shape={kind='all'}, reduction={'average','average'},
    time_control={min={iter=0}, max={iter=maximum_iterations}, interval={iter=1}},
    output={format='ascii'}},
  {label='pdf_safety', folder='tracking/', variable={'pdf'},
    shape={kind='all'}, reduction={'min'},
    time_control={min={iter=0}, max={iter=maximum_iterations}, interval={iter=1}},
    output={format='ascii'}}
}
