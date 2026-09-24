-- Fully periodic 4^3 force-only binary oracle.
folder = 'mesh/'
timing_file = 'seeder_timing.res'
comment = 'force-only oracle; no wall/inlet/outlet/pressure/adaptive boundary'
logging = { level = 3 }
minlevel = 2
bounding_cube = { origin = {0.0, 0.0, 0.0}, length = 4.6746357867857693e-07 }
e = 2.7863000552569445e-14
s = 4.6746363440457801e-07
L = 4.6746357867857693e-07
spatial_object = {
  { attribute = {kind='seed'}, geometry = {kind='canoND', object={origin={L/2,L/2,L/2}}} },
  { attribute = {kind='periodic'}, geometry = {kind='periodic', object={
    plane1={origin={-e,-e,-e}, vec={{0,s,0},{0,0,s}}},
    plane2={origin={L+e,-e,-e}, vec={{0,0,s},{0,s,0}}}
  }} },
  { attribute = {kind='periodic'}, geometry = {kind='periodic', object={
    plane1={origin={-e,-e,-e}, vec={{0,0,s},{s,0,0}}},
    plane2={origin={-e,L+e,-e}, vec={{s,0,0},{0,0,s}}}
  }} },
  { attribute = {kind='periodic'}, geometry = {kind='periodic', object={
    plane1={origin={-e,-e,-e}, vec={{s,0,0},{0,s,0}}},
    plane2={origin={-e,-e,L+e}, vec={{0,s,0},{s,0,0}}}
  }} }
}
