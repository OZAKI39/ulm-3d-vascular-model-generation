-- Generated production configuration; official Seeder syntax.
folder = 'mesh/'
comment = 'ROI003274 uniform 0.20 um CFD lattice'
debug = { debugMode = false, debugFiles = false, debugMesh = 'debug/' }
minlevel = 9
bounding_cube = { origin = { 7.8082771579034594e-05, 5.4076359424537688e-05, 8.8206398940868488e-05 }, length = 0.0001024 }
spatial_object = {
  {
    attribute = { kind = 'boundary', label = 'wall', level = minlevel, calc_dist = true },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/wall.stl' }}
  },
  {
    attribute = { kind = 'boundary', label = 'inlet', level = minlevel },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/inlet.stl' }}
  },
  {
    attribute = { kind = 'boundary', label = 'outlet_01', level = minlevel },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/outlet_01.stl' }}
  },
  {
    attribute = { kind = 'boundary', label = 'outlet_02', level = minlevel },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/outlet_02.stl' }}
  },
  {
    attribute = { kind = 'boundary', label = 'outlet_03', level = minlevel },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/outlet_03.stl' }}
  },
  {
    attribute = { kind = 'seed' },
    geometry = { kind = 'canoND', object = { origin = { 0.00010367370500008977, 5.75542332387685e-05, 0.00015801755415184839 } } }
  }
}
