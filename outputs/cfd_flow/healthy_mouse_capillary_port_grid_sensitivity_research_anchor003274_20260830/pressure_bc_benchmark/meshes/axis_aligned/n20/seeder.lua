-- Isolated pressure-BC pipe benchmark; not production CFD.
folder = 'mesh/'
comment = 'fixed physical Hagen-Poiseuille pipe, N=20'
debug = { debugMode = false, debugFiles = false, debugMesh = 'debug/' }
minlevel = 8
bounding_cube = { origin = { -2.0194426598914522e-05, -2.0194426598914522e-05, -2.0194426598914522e-05 }, length = 4.0388853197829045e-05 }
spatial_object = {
  { attribute = { kind = 'boundary', label = 'wall', level = minlevel, calc_dist = true }, geometry = { kind = 'stl', object = { filename = 'geometry/wall.stl' } } },
  { attribute = { kind = 'boundary', label = 'inlet', level = minlevel }, geometry = { kind = 'canoND', object = { origin = { 1.735458535844217e-06, -1.735458535844217e-06, -9.4661374682411823e-06 }, vec = { { -3.4709170716884339e-06, 0, 0 }, { 0, 3.4709170716884339e-06, 0 } } } } },
  { attribute = { kind = 'boundary', label = 'outlet', level = minlevel }, geometry = { kind = 'canoND', object = { origin = { 1.735458535844217e-06, -1.735458535844217e-06, 9.4661374682411823e-06 }, vec = { { 0, 3.4709170716884339e-06, 0 }, { -3.4709170716884339e-06, 0, 0 } } } } },
  { attribute = { kind = 'seed' }, geometry = { kind = 'canoND', object = { origin = { 0, 0, 0 } } } }
}
