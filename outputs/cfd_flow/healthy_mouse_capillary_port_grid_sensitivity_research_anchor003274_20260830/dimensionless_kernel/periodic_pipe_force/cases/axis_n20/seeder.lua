-- Frozen dimensionless-kernel periodic Pipe Force mesh.
folder = 'mesh/'
timing_file = 'seeder_timing.res'
comment = 'axis_n20; wall_libb continuous q; periodic axial direction'
logging = { level = 3 }
minlevel = 8
bounding_cube = { origin = { 0.0, 0.0, 0.0 }, length = 4.0388853197829045e-05 }
spatial_object = {
  { attribute = { kind = 'seed' },
     geometry = { kind = 'canoND', object = { origin = { 2.0194426598914522e-05, 2.0194426598914522e-05, 2.0194426598914522e-05 } } } },
  { attribute = { kind = 'boundary', label = 'wall', level = minlevel,
                     calc_dist = true, flood_diagonal = false },
     geometry = { kind = 'stl', object = { filename = 'geometry/wall.stl' } } },
  { attribute = { kind = 'periodic' },
     geometry = { kind = 'periodic', object = {
       plane1 = { origin = { 2.208765409256276e-05, 1.8301199105266284e-05, 1.0728250612861377e-05 },
                  vec = { { -3.7864549872964729e-06, 0, 0 }, { 0, 3.7864549872964729e-06, 0 } } },
       plane2 = { origin = { 2.208765409256276e-05, 1.8301199105266284e-05, 2.9660602584967669e-05 },
                  vec = { { 0, 3.7864549872964729e-06, 0 }, { -3.7864549872964729e-06, 0, 0 } } }
     } } }
}
