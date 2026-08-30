-- Axis-aligned ideal numerical inlet plane preflight; Seeder only.
folder = 'mesh/'
comment = 'ROI003274 ideal numerical inlet plane, uniform 0.20 um lattice'
debug = { debugMode = false, debugFiles = false, debugMesh = 'debug/' }
minlevel = 9
bounding_cube = { origin = { 78.082771579034599, 54.076359424537685, 88.206398940868482 }, length = 102.39999999999999 }
spatial_object = {
  {
    attribute = { kind = 'boundary', label = 'wall', level = minlevel, calc_dist = true },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/wall.stl' } }
  },
  {
    attribute = { kind = 'boundary', label = 'inlet', level = minlevel },
    geometry = { kind = 'canoND', object = { origin = { 101.68277157903459, 55.476359424537691, 158.80639894086849 }, vec = { { 3.9999999999999889, 0, 0 }, { 0, 4.2000000000000002, 0 } } } }
  },
  {
    attribute = { kind = 'boundary', label = 'outlet_01', level = minlevel },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/outlet_01.stl' } }
  },
  {
    attribute = { kind = 'boundary', label = 'outlet_02', level = minlevel },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/outlet_02.stl' } }
  },
  {
    attribute = { kind = 'boundary', label = 'outlet_03', level = minlevel },
    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/outlet_03.stl' } }
  },
  {
    attribute = { kind = 'seed' },
    geometry = { kind = 'canoND', object = { origin = { 103.67370500008977, 57.554233238768504, 158.0175541518484 } } }
  }
}
