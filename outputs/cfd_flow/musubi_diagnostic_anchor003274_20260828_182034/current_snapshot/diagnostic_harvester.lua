-- Diagnostic harvest of the non-converged iteration-162464 restart.
require 'musubi'
restart = { read = 'restart/roi003274_steady_lbm_lastHeader.lua' }
tracking = {
  label = 'current_snapshot_162464',
  folder = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/musubi_diagnostic_anchor003274_20260828_182034/current_snapshot/raw_vtk/',
  variable = { 'pressure_phy', 'velocity_phy', 'vel_mag_phy' },
  shape = { kind = 'all' },
  output = { format = 'vtk' }
}
