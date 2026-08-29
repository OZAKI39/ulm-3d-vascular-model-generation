-- Isolated tiny Cartesian z-channel for adaptive flux validation.
folder = 'mesh/'
logging = { level = 4 }
debug = { debugMode = false, debugFiles = false }
dx = 2.0e-7
width = 8*dx
length = 16*dx
level = 5
cube = (2^level)*dx
eps = cube/(2^22)
bounding_cube = { origin = {-dx, -dx, -dx}, length = cube }
minlevel = level
spatial_object = {
  { attribute = {kind='seed'}, geometry = {kind='canoND', object={origin={width/2,width/2,length/2}}}},
  { attribute = {kind='boundary',label='inlet'}, geometry={kind='canoND',object={origin={-eps,-eps,-eps},vec={{width+2*eps,0,0},{0,width+2*eps,0}}}}},
  { attribute = {kind='boundary',label='outlet'}, geometry={kind='canoND',object={origin={-eps,-eps,length+eps},vec={{width+2*eps,0,0},{0,width+2*eps,0}}}}},
  { attribute = {kind='boundary',label='wall_x0'}, geometry={kind='canoND',object={origin={-eps,-eps,-eps},vec={{0,width+2*eps,0},{0,0,length+2*eps}}}}},
  { attribute = {kind='boundary',label='wall_x1'}, geometry={kind='canoND',object={origin={width+eps,-eps,-eps},vec={{0,width+2*eps,0},{0,0,length+2*eps}}}}},
  { attribute = {kind='boundary',label='wall_y0'}, geometry={kind='canoND',object={origin={-eps,-eps,-eps},vec={{width+2*eps,0,0},{0,0,length+2*eps}}}}},
  { attribute = {kind='boundary',label='wall_y1'}, geometry={kind='canoND',object={origin={-eps,width+eps,-eps},vec={{width+2*eps,0,0},{0,0,length+2*eps}}}}}
}
