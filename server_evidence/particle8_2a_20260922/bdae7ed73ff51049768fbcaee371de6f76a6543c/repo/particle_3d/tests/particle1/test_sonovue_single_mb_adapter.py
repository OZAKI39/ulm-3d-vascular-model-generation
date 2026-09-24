import numpy as np
import pytest
from particle_3d.sonovue_adapter import read_sonovue,sample_single_validation_size
from particle_3d.particle1_audit import VALIDATION_SEED


def test_exact_original_sampler_call_and_reproducibility(sonovue_root,size):
    contract,distribution,_=read_sonovue(sonovue_root)
    direct=distribution.sample_diameters(n=1,seed=VALIDATION_SEED)
    assert size.diameter_um == direct[0]
    assert sample_single_validation_size(sonovue_root,seed=VALIDATION_SEED)==size
    assert size.diameter_m==size.diameter_um*1e-6 and size.radius_m==.5*size.diameter_m
    assert size.histogram_sha256==contract['histogram_sha256']
    assert size.sampler_version==contract['sampler_version']
    assert size.N==1 and size.formal_simulation_population is False
    assert size.role=='DEMO / VALIDATION ONLY'
    assert size.numpy_version==np.__version__


def test_seed_is_explicit_and_invalid_seed_rejected(sonovue_root):
    with pytest.raises(TypeError): sample_single_validation_size(sonovue_root)
    with pytest.raises(ValueError): sample_single_validation_size(sonovue_root,seed=-1)
