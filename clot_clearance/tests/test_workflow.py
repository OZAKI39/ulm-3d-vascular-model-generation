import copy,json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.runner import simulate


@pytest.fixture
def config():
    root=Path(__file__).resolve().parents[1]
    c=json.loads((root/'configs/straight_pipe.json').read_text())
    c['clot']['cells']=[6,4,4];c['clot']['origin_m']=[-.000375,-.00025,-.0011]
    c['streaming']['bubble_center_m']=[0,0,-.0005]
    c['simulation']['number_of_macro_steps']=3
    c['simulation']['representative_frequency_Hz']=1000.
    c['simulation']['warmup_cycles']=0
    c['streaming']['include_pipe_traction']=False
    return c


def test_A_no_streaming(config):
    config['streaming']['traction_scale_Pa']=0
    r=simulate(config,quiet=True)
    assert r['summary']['final']['maximum_displacement_um']<1e-7
    assert r['integrity'].min()==1
    assert r['summary']['final']['number_of_fragments']==1


def test_B_small_uniform_load(config):
    config['streaming'].update(traction_scale_Pa=1.,uniform=True)
    r=simulate(config,quiet=True)
    assert r['summary']['final']['maximum_displacement_um']>0.001
    assert r['integrity'].min()==1
    assert r['summary']['final']['fixed_base_error_m']==0


def test_C_D_E_local_damage_monotonic_cycles_and_zero_rate(config):
    config['streaming']['traction_scale_Pa']=160.
    config['damage'].update(Q_ref=.004,C_damage=.000002)
    r=simulate(config,quiet=True);h=r['history']
    assert h[-1]['mean_particle_damage']>0
    assert np.all(np.diff([z['mean_particle_damage'] for z in h])>=-1e-15)
    # Localized field damage is greater in upper loaded layer than fixed base.
    g=r['cloud'];D=r['states'][-1]['damage']
    assert D[g.X[:,2]>g.X[:,2].mean()].mean()>D[g.fixed].mean()
    config['damage']['C_damage']=0
    z=simulate(config,quiet=True)
    assert np.all(z['integrity']==1)
    config['damage']['enabled']=False
    z=simulate(config,quiet=True)
    assert np.all(z['integrity']==1)


def test_safety_bad_dt_and_preserve_existing_output(config,tmp_path):
    config['simulation']['dt_s']=.1
    with pytest.raises(ValueError,match='Explicit dt'):simulate(config,quiet=True)
    config['simulation']['dt_s']=.000002
    with pytest.raises(FileExistsError):simulate(config,out=tmp_path,quiet=True)
