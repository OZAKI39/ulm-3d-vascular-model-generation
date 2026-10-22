import numpy as np,pytest
from scripts.straight_tube_exp1.f10z.prepare import generate,action
from scripts.straight_tube_exp1.f10z.state import series,validate_trajectory
from mevo_mb.control.scheduled_pzt import Schedule,Segment
from scipy.spatial.distance import pdist

def test_initial_legality_and_determinism():
 a,_=generate();b,_=generate()
 for k,x in a.items():
  assert np.array_equal(x,b[k]);assert pdist(x).min()>=3.56e-6;assert np.linalg.norm(x[:,1:],axis=1).max()+1.78e-6<.0008
  assert np.abs(x).max()<54e-6
@pytest.mark.parametrize('name',['OFF','XM','XP','PAIR','YM','YP','AMP-5','PHASE+3','PAIR_OFF'])
def test_schedules(name):
 s=action(name);assert Schedule.from_record(s.record())==s
 steps=list(s.steps(.7e-6));assert abs(sum(v[1] for v in steps)-s.duration_s)<1e-13
 assert all(v[2].frequency_Hz==2.5e6 for v in steps)
 if name=='PAIR_OFF':assert any(abs(v[0]-.001)<1e-14 and v[2].segment_id==1 for v in steps)
def test_range_rejects():
 with pytest.raises(ValueError):Segment(0,0,.001,'ONE_HOT_PZT',(1501,0,0,0),(0,)*4,(True,False,False,False))
def test_translation_not_convergence():
 x=generate()[0]['medium_train'];q=np.array([x,x+[1e-7,0,0]]);m=series(np.array([0.,.001]),q)
 assert abs(m['rg_m'][1]-m['rg_m'][0])<1e-20
 assert np.max(abs(m['relative_velocity_m_s']))<1e-16
 assert abs(m['radial_mean_m_s'][0])<1e-16
@pytest.mark.parametrize('bad',['time','ids','shape','nan'])
def test_bad_trajectory(bad):
 x=np.zeros((2,32,3));t=np.array([0.,.001]);ids=np.arange(32)
 if bad=='time':t[:]=0
 if bad=='ids':ids=ids[::-1]
 if bad=='shape':x=x[:,:31]
 if bad=='nan':x[0,0,0]=np.nan
 with pytest.raises(ValueError):validate_trajectory(t,x,ids)
def test_stationary_and_incomplete():
 x=generate()[0]['medium_train'];m=series(np.array([0.]),x[None]);assert len(m['velocity_m_s'])==0 and len(m['time_s'])==1
 m=series(np.array([0.,.0007,.001]),np.repeat(x[None],3,0));assert np.max(abs(m['velocity_m_s']))==0

def test_f10y_distribution_definitions():
 import importlib.util
 from scripts.straight_tube_exp1.f10z.common import ROOT
 p=ROOT.parents[1]/'external_connect/vision_tracking_research/realtime_collective_observation_v1/src/aggregation.py'
 spec=importlib.util.spec_from_file_location('frozen_f10y_distribution',p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 x=generate()[0]['medium_train'];reference=module.distribution(x[:,:2]*1e6);m=series(np.array([0.]),x[None])
 np.testing.assert_allclose(reference['center'],x[:,:2].mean(0)*1e6,atol=1e-13)
 np.testing.assert_allclose(reference['rg_px'],m['xy_rg_m'][0]*1e6)
 np.testing.assert_allclose(reference['shape_tensor'],m['xy_tensor_m2'][0]*1e12)
 np.testing.assert_allclose(reference['nn_px'],m['xy_nn_m'][0]*1e6)
 # Coordinates were scaled to micrometres for formula comparison only, not camera pixels.

def test_failed_trajectory_export_does_not_fill_time(tmp_path):
 import h5py,json
 from scripts.straight_tube_exp1.f10z.state import export
 x=generate()[0]['medium_train'];folder=tmp_path/'failed';folder.mkdir()
 with h5py.File(folder/'trajectory.h5','w') as h:
  h['time_s']=[0.,.0007];h['positions_m']=np.array([x,x]);h['particle_id']=np.arange(32);h.attrs['schedule_json']=json.dumps(action('PAIR').record())
 (folder/'RESULT.json').write_text(json.dumps({'requested_s':.002,'completed_s':.0007,'status':'FAIL'}))
 out=tmp_path/'state.h5';export(folder,out)
 with h5py.File(out) as h:
  assert len(h['time_s'])==2 and h['time_s'][-1]==.0007
  assert h.attrs['requested_s']==.002 and h.attrs['completed_s']==.0007 and h.attrs['no_padding']
