-- Arbitrary-q wall_libb oracle: q=0.25 and q=0.75 in one mesh.
folder = 'mesh/'
timing_file = 'seeder_timing.res'
comment = 'wall_libb arbitrary-q binary oracle; periodic y/z'
logging = {level=3}
minlevel = 2
bounding_cube = {origin={0.0,0.0,0.0}, length=4.6746357867857693e-07}
e=2.7863000552569445e-14; s=4.6746363440457801e-07; L=4.6746357867857693e-07
spatial_object = {
  {attribute={kind='seed'}, geometry={kind='canoND', object={origin={2*L/4,L/2,L/2}}}},
  {attribute={kind='boundary', label='wall', level=minlevel, calc_dist=true, flood_diagonal=false},
    geometry={kind='stl', object={filename='geometry/wall.stl'}}},
  {attribute={kind='periodic'}, geometry={kind='periodic', object={
    plane1={origin={-e,-e,-e}, vec={{0,0,s},{s,0,0}}},
    plane2={origin={-e,L+e,-e}, vec={{s,0,0},{0,0,s}}}
  }}},
  {attribute={kind='periodic'}, geometry={kind='periodic', object={
    plane1={origin={-e,-e,-e}, vec={{s,0,0},{0,s,0}}},
    plane2={origin={-e,-e,L+e}, vec={{0,s,0},{s,0,0}}}
  }}}
}
