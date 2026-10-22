import sys,json,copy,time
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from xcommon import parameters,fit_affine
from fit_adapter import fit_affine as fast
from channel import Channel
from observer import groups_from_candidates,finite_list
from channel_state import validate_channel,to_v1,to_v2,expire
from scipy.linalg import expm
def test_six_coordinate_reference():
 A=np.array([[.1,-.2],[.3,-.05]]);center=np.array([750.,410.]);velocity=np.array([15.,-8.]);L=np.zeros((3,3));L[:2,:2]=A;L[:2,2]=velocity-A@center;H=expm(L*.0999)
 p=parameters(H[:2,:2],H[:2,2],.0999,center);np.testing.assert_allclose(p,[15,-8,.05,.5,.075,.05],atol=1e-10)
 shift=np.array([650.,300.]);b=H[:2,2]+(H[:2,:2]-np.eye(2))@shift
 np.testing.assert_allclose(parameters(H[:2,:2],b,.0999,center-shift),p,atol=1e-10)
def test_window_coverage_not_center_only():
 c=Channel();assert np.all(c.coverage[c.safe]>=1.)
 center=c.local[c.y.astype(int),c.x.astype(int)]>0
 assert np.count_nonzero(center&~c.safe)>0
def test_nulls_and_no_fake_zeros():
 assert finite_list([np.nan,1.,np.inf])==[None,1.,None]
def test_fast_fit_frozen_gates():
 rng=np.random.default_rng(324);x=rng.uniform(25,215,(60,2));d=np.c_[.02*x[:,0]+2,-.01*x[:,1]-1]+rng.normal(0,.2,(60,2))
 a=fit_affine(x,d,.1,(240,240),'legacy');b=fast(x,d,.1,(240,240),'fast')
 np.testing.assert_allclose(a['raw_estimate'],b['raw_estimate'],atol=1e-5);assert a['valid']==b['valid']
def test_degenerate_support_rejected():
 x=np.c_[np.arange(10),np.ones(10)];r=fast(x,np.ones_like(x),.1,(240,240),'test');assert not any(r['valid'].values());assert all(v is None for v in r['values'].values())
def test_state_legacy_adapters_and_expiry():
 from xcommon import empty
 state=dict(schema_version='M4-channel-3.0',source='FILE',session_id='TEST',frame_id=4,source_PTS=.1332,host_receive_timestamp=100.,ROI=[432,0,912,1120],calibration_status='PIXEL_ONLY_NO_PHYSICAL_SCALE',stale_or_discontinuous=False,channel_mask_valid=True,backend='B4',training_ready=False,hardware_control_executed=False,time_window={},groups=[],windows=dict(x=[700],y=[400],vx=[None],vy=[None],valid=[False],fresh=[False],stale=[True],safe=[True],last_update_pts=[None],age_ms=[None]))
 assert validate_channel(state)
 v1=to_v1(state);v2=to_v2(state);assert not any(v1['motion']['valid'].values());assert v1['schema_version']=='M4-collective-1.0'
 import jsonschema
 schema=json.loads((ROOT.parent/'collective_motion_v2/configs/collective_state_v2.schema.json').read_text());jsonschema.validate(v2,schema)
 stale=expire(state,100.101);assert stale['stale_or_discontinuous'];assert stale['windows']['vx']==[None]

def test_robust_jacobian_is_not_reused_after_scipy_mutates_it():
    import numpy as np
    from xcommon import fit_affine as old
    from fit_adapter import fit_affine as new
    rng=np.random.default_rng(930)
    xy=rng.uniform(5,235,(50,2))
    displacement=np.tile([2.,-1.],(50,1))+rng.normal(0,.05,(50,2))
    displacement[:9]+=rng.normal(0,3,(9,2))
    a=old(xy,displacement,.0999,(240,240),'old')
    b=new(xy,displacement,.0999,(240,240),'new')
    assert a['valid']==b['valid']
    np.testing.assert_allclose(a['raw_estimate'],b['raw_estimate'],atol=1e-4,rtol=0)
