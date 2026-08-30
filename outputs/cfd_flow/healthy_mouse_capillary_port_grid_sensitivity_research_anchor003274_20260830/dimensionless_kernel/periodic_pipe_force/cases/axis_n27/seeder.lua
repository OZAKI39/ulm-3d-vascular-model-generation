-- Frozen dimensionless-kernel periodic Pipe Force mesh.
folder = 'mesh/'
timing_file = 'seeder_timing.res'
comment = 'axis_n27; wall_libb continuous q; periodic axial direction'
logging = { level = 3 }
minlevel = 8
bounding_cube = { origin = { 0.0, 0.0, 0.0 }, length = 2.9917669035428923e-05 }
spatial_object = {
  { attribute = { kind = 'seed' },
     geometry = { kind = 'canoND', object = { origin = { 1.4958834517714462e-05, 1.4958834517714462e-05, 1.4958834517714462e-05 } } } },
  { attribute = { kind = 'boundary', label = 'wall', level = minlevel,
                     calc_dist = true, flood_diagonal = false },
     geometry = { kind = 'stl', object = { filename = 'geometry/wall.stl' } } },
  { attribute = { kind = 'periodic' },
     geometry = { kind = 'periodic', object = {
       plane1 = { origin = { 1.6770255885093947e-05, 1.3147413150334976e-05, 5.4926685177607135e-06 },
                  vec = { { -3.6228427347589708e-06, 0, 0 }, { 0, 3.6228427347589708e-06, 0 } } },
       plane2 = { origin = { 1.6770255885093947e-05, 1.3147413150334976e-05, 2.4425000517668212e-05 },
                  vec = { { 0, 3.6228427347589708e-06, 0 }, { -3.6228427347589708e-06, 0, 0 } } }
     } } }
}
