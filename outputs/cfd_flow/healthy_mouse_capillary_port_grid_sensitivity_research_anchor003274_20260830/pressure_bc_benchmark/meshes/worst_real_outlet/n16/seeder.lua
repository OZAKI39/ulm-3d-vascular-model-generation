-- Isolated pressure-BC pipe benchmark; not production CFD.
folder = 'mesh/'
comment = 'fixed physical Hagen-Poiseuille pipe, N=16'
debug = { debugMode = false, debugFiles = false, debugMesh = 'debug/' }
minlevel = 8
bounding_cube = { origin = { -2.5243033248643153e-05, -2.5243033248643153e-05, -2.5243033248643153e-05 }, length = 5.0486066497286306e-05 }
spatial_object = {
  { attribute = { kind = 'boundary', label = 'wall', level = minlevel, calc_dist = true }, geometry = { kind = 'stl', object = { filename = 'geometry/wall.stl' } } },
  { attribute = { kind = 'boundary', label = 'inlet', level = minlevel }, geometry = { kind = 'canoND', object = { origin = { 1.561097935389852e-06, -8.6724324617025868e-06, -4.2406697183625915e-06 }, vec = { { -3.4703302055271554e-06, 5.1362625181231359e-08, 3.7887516130751122e-08 }, { 0, 2.0604025528482922e-06, -2.7932072316164535e-06 } } } } },
  { attribute = { kind = 'boundary', label = 'outlet', level = minlevel }, geometry = { kind = 'canoND', object = { origin = { 1.9092322701373034e-06, 6.5606672836730623e-06, 6.9959894338482939e-06 }, vec = { { 0, 2.0604025528482922e-06, -2.7932072316164535e-06 }, { -3.4703302055271554e-06, 5.1362625181231359e-08, 3.7887516130751122e-08 } } } } },
  { attribute = { kind = 'seed' }, geometry = { kind = 'canoND', object = { origin = { 0, 0, 0 } } } }
}
