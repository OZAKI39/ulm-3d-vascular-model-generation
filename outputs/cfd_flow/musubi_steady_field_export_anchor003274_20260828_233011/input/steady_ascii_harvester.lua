-- One-shot full field export from the frozen project-steady restart.
package.path = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/musubi_project_steady_confirmation_anchor003274_20260828_225334/?.lua;' .. package.path
require 'diagnostic_musubi'

restart = {
  read = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/musubi_project_steady_confirmation_anchor003274_20260828_225334/restart/roi003274_steady_lbm_lastHeader.lua'
}

tracking = {
  label = 'f',
  folder = '/tmp/u3d/x/',
  variable = { 'pressure_phy', 'velocity_phy' },
  shape = { kind = 'all' },
  output = { format = 'asciiSpatial', use_get_point = false }
}
