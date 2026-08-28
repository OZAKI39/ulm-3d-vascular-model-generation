-- One final full-volume harvest from the terminal restart.
require 'musubi'
restart = { read = 'restart/roi003274_steady_lbm_lastHeader.lua' }
tracking = {
  label = 'flow_field',
  folder = '../flow/',
  variable = { 'pressure_phy', 'velocity_phy' },
  shape = { kind = 'all' },
  output = { format = 'vtk' }
}
