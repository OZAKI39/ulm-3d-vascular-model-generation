import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import copy
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from piv_display import PIVDisplayFilter,piv_arrows,display_vectors,display_is_smoothed
from pipeline import observation_record


def state(t,velocity=20.,x=24.,session='s',revision=0):
    return dict(valid=True,arrival_monotonic=t,session_id=session,geometry_revision=revision,
        vectors=np.array([[x,24.,velocity,0.]]),geometry=None)


def test_fixed_anchors_reduce_display_jitter_without_changing_measurements():
    f=PIVDisplayFilter();raw=[];shown=[];positions=[]
    for i in range(40):
        s=state(1.+i*.05,20.+(-1)**i*6,x=24.+(i%2)*24)
        original=s['vectors'].copy();v=f.update(s)
        assert np.array_equal(original,s['vectors'])
        raw.append(original[0,2]);shown.append(v[0,2]);positions.append(v[0,:2])
    assert np.max(np.ptp(positions,axis=0))==0
    assert np.std(np.diff(shown[5:]))<.4*np.std(np.diff(raw[5:]))
    assert np.allclose(np.mean(shown[5:]),20,atol=.3)


def test_missing_support_invalidity_gap_and_context_reset_no_ghost_arrows():
    f=PIVDisplayFilter();f.update(state(1.))
    empty=state(1.05);empty['vectors']=np.empty((0,4))
    assert not len(f.update(empty))
    assert f.update(state(1.1,50.))[0,2]==50.
    invalid=state(1.15);invalid['valid']=False
    assert not len(f.update(invalid))
    assert f.update(state(1.2,10.))[0,2]==10.
    assert f.update(state(1.7,50.))[0,2]==50. # Original TTL exceeded.
    assert f.update(state(1.75,10.,session='new'))[0,2]==10.
    assert f.update(state(1.8,50.,session='new',revision=2))[0,2]==50.
    assert f.update(state(1.6,10.,session='new',revision=2))[0,2]==10. # Clock reset.


def test_stop_and_reversal_are_not_delayed_by_filter():
    f=PIVDisplayFilter();f.update(state(1.,30.))
    stopped=f.update(state(1.05,0.));assert not piv_arrows(stopped,sample=False)
    f.update(state(1.1,30.))
    reverse=f.update(state(1.15,-20.))
    assert reverse[0,2]==-20.


def test_coarse_field_stays_inside_channel_and_exclusions():
    s=state(1.);s['vectors']=np.array([[48,48,20,0],[80,48,20,0]])
    s['geometry']=dict(image_wh=[144,144],crop_xywh=[0,0,144,144],
        channel_polygon=[[0,0],[70,0],[70,100],[0,100]],excluded_polygons=[])
    v=PIVDisplayFilter().update(s)
    assert len(v)==1 and np.array_equal(v[0,:2],[36,36])
    s['geometry']['excluded_polygons']=[[[30,30],[40,30],[40,40],[30,40]]]
    assert not len(PIVDisplayFilter().update(s))


def test_raw_inspection_and_scientific_record_remain_unfiltered():
    s=state(1.);s['piv_display_vectors']=PIVDisplayFilter().update(s)
    original=copy.deepcopy(s)
    assert display_is_smoothed(s,{'Smooth PIV':True})
    assert not display_is_smoothed(s,{'Smooth PIV':False})
    assert np.array_equal(display_vectors(s,False),s['vectors'])
    assert np.array_equal(display_vectors(s,True),s['piv_display_vectors'])
    record=observation_record(s)
    assert 'piv_display_vectors' not in record and np.array_equal(record['vectors'],original['vectors'])
