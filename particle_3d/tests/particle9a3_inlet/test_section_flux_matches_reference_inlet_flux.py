import numpy as np
from particle_3d.interior_section import evaluate_candidate, cut_tetrahedra, p1_flux


def test_constant_velocity_analytic_area_flux(tube):
    grid, candidate, topo, wall, caps = tube
    result, _ = evaluate_candidate(grid, candidate, topo, wall, caps, 2e-15, 1e-6)
    assert result['accepted']
    assert abs(result['area_m2']/1e-12-1) < 1e-13
    assert result['flux_relative_error'] < 1e-13


def test_positive_part_linear_fem_integrated_exactly(tube):
    grid, candidate, *_ = tube
    # u_z = 4000*x; positive half-square integral = 4000*(0.5e-6)^2/2*1e-6.
    grid['Velocity'][:, 2] = 4000*grid.points[:, 0]
    section = cut_tetrahedra(grid, candidate['center_m'], candidate['normal'])
    flux = p1_flux(section, candidate['normal'])
    assert abs(flux['signed_Q_m3_s']) < 1e-28
    assert abs(flux['positive_Q_m3_s']/5e-16-1) < 1e-13


def test_percent_flux_deficit_cannot_be_accepted(tube):
    grid, candidate, topo, wall, caps = tube
    grid['Velocity'] *= .99
    result, _ = evaluate_candidate(grid, candidate, topo, wall, caps, 2e-15, 1e-6)
    assert result['geometry_pass']
    assert result['rejected_reason'] == 'SECTION_FLUX_MISMATCH'
    assert not result['accepted']
