import copy
import json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.geometry import make_cloud
from pd_clot.regularization.damage import EnergyBondDamage, LegacyDamageWithWarning


def fixture():
    root=Path(__file__).parents[1]
    c=json.loads((root/'configs/streaming_regularized_demo_coarse.json').read_text())
    cal=json.loads((root/'verification/regularization/calibration_coarse/CALIBRATION.json').read_text())
    return c,cal,make_cloud(c['clot'])


def test_calibration_identity_and_full_damage_required():
    c,cal,g=fixture();d=EnergyBondDamage(g,c,cal)
    d.advance(np.full(len(g.pairs),1e-6),1000)
    assert d.active.all() and np.all(d.D>0)
    d.advance(np.ones(len(g.pairs)),1000)
    assert not d.active.any() and np.all(d.D==1)
    d.advance(np.zeros(len(g.pairs)),1000)
    assert np.all(d.D==1) and np.all(d.g==0)
    wrong=copy.deepcopy(cal);wrong['spacing_m']*=.5
    with pytest.raises(ValueError,match='mismatched'):EnergyBondDamage(g,c,wrong)
    c['regularized_damage']['D_break']=.25
    with pytest.raises(ValueError,match='D_break'):EnergyBondDamage(g,c,cal)


def test_legacy_warning_retains_old_damage_law():
    c,_,_=fixture();cfg={'mode':'cyclic_accumulation','driver':'bond_stretch_amplitude','D_break':.25,'DeltaN':1000,'C_damage':.001,'Q_ref':.01,'m_damage':2,'enabled':True}
    with pytest.warns(RuntimeWarning,match='forced_failure_verification'):
        d=LegacyDamageWithWarning(2,cfg)
    d.advance(np.array([0.,.02]));assert d.active.tolist()==[True,False]
