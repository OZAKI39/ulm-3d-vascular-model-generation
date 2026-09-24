-- Frozen dimensionless-kernel periodic Pipe Force mesh.
folder = 'mesh/'
timing_file = 'seeder_timing.res'
comment = 'axis_n16; wall_libb continuous q; periodic axial direction'
logging = { level = 3 }
minlevel = 8
bounding_cube = { origin = { 0.0, 0.0, 0.0 }, length = 5.0486066497286306e-05 }
spatial_object = {
  { attribute = { kind = 'seed' },
     geometry = { kind = 'canoND', object = { origin = { 2.5243033248643153e-05, 2.5243033248643153e-05, 2.5243033248643153e-05 } } } },
  { attribute = { kind = 'boundary', label = 'wall', level = minlevel,
                     calc_dist = true, flood_diagonal = false },
     geometry = { kind = 'stl', object = { filename = 'geometry/wall.stl' } } },
  { attribute = { kind = 'periodic' },
     geometry = { kind = 'periodic', object = {
       plane1 = { origin = { 2.7215145221193397e-05, 2.3270921276092908e-05, 1.5776847633137012e-05 },
                  vec = { { -3.9442239451004922e-06, 0, 0 }, { 0, 3.9442239451004922e-06, 0 } } },
       plane2 = { origin = { 2.7215145221193397e-05, 2.3270921276092908e-05, 3.4709218864149293e-05 },
                  vec = { { 0, 3.9442239451004922e-06, 0 }, { -3.9442239451004922e-06, 0, 0 } } }
     } } }
}
