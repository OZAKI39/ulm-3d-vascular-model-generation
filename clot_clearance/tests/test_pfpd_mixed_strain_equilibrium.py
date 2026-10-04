"""Evidence gates for the requested twelve frozen-phase equilibrium experiments."""
import numpy as np
from scripts.spectral_split_audit import OUT,read,ca

def test_all_requested_intact_equilibria_pass_original_force_tolerance():
    rows=read(OUT/'equilibrium/ALL.json')
    assert {(r['level'],r['strain']) for r in rows}=={(f'C{n}',e) for n in [3,5,7] for e in [.00025,.005,.009,.01]}
    assert all(r['accepted'] and all(r['gates'].values()) and r['force_residual']<=1e-9 for r in rows)

def test_intact_experiment_never_updates_phase():
    for r in read(OUT/'equilibrium/ALL.json'):
        arrays=ca.read_npz(OUT/f'equilibrium/{r["level"]}/{r["strain"]:.8f}/STATE.npz')
        assert np.count_nonzero(arrays['s'])==0 and r['gates']['phase_frozen']

def test_fine_equilibrium_is_mixed_with_positive_J():
    for r in read(OUT/'equilibrium/ALL.json'):
        assert r['minimum_J']>0
        if r['level'] in ['C5','C7']:
            assert r['min_principal_green_strain']<-1e-9
            assert r['negative_principal_bond_fraction']>0 and r['Psi_minus_fraction']>0

def test_negative_strain_fraction_is_not_compression_energy_fraction():
    for r in read(OUT/'equilibrium/ALL.json'):
        np.testing.assert_allclose(r['total_Psi_plus_J']+r['total_Psi_minus_J'],r['elastic_energy_J'],rtol=2e-12)
        np.testing.assert_allclose(r['Psi_plus_fraction']+r['Psi_minus_fraction'],1,rtol=2e-15)
        if r['level'] in ['C5','C7']:assert r['negative_principal_bond_fraction']>10*r['Psi_minus_fraction']

def test_compressive_locations_stay_inside_fixed_physical_domain():
    for r in read(OUT/'equilibrium/ALL.json'):
        for key in ['most_compressive_particle','most_compressive_BAQP']:
            point=r[key];assert point['clamp_distance_m']>=0 and point['side_distance_m']>=0
            assert np.max(abs(np.array(point['position_m'])))<=.000125
