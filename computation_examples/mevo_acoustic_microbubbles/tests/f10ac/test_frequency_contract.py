import json,numpy as np,h5py,pytest
from mevo_mb.control.scheduled_pzt import Segment,Schedule
from mevo_mb.fast_sim.nearfield import F2Physics
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import FrequencySegment,FixedResponse,make_schedule

def test_legacy_frequency_remains_rejected():
 with pytest.raises(ValueError):Segment(0,0,.001,'FOUR_CHANNEL_INDEPENDENT',(0,)*4,(0,)*4,(False,)*4,frequency_Hz=2400000.)

def test_frequency_step_clock_and_identity():
 s=make_schedule([0,0,500,500],[0,0,0,1],.0000025,2400000.)
 rows=list(s.steps(1e-6));assert np.allclose([r[0] for r in rows],[0,1e-6,2e-6],atol=1e-16)
 assert sum(r[1] for r in rows)==pytest.approx(.0000025)
 assert all(r[2].event()['frequency_Hz']==2400000. for r in rows)

@pytest.mark.parametrize('freq',[0,2499999.,2800000.,float('nan')])
def test_no_frequency_extrapolation(freq):
 with pytest.raises(ValueError):make_schedule([100]*4,[0]*4,frequency=freq)

def test_strict_pressure_and_phase(tmp_path):
 path=tmp_path/'r.h5';v=np.tile(np.arange(1,9),(3,1))*(1+2j)*np.array([1,2,3])[:,None]
 with h5py.File(path,'w') as h:
  h['pressure_Pa']=[1,2,3];h['harmonics_m3']=v;h.attrs['metadata_json']=json.dumps(dict(frequency_Hz=2400000.,radius_m=1.78e-6,validated=True))
 r=FixedResponse(path,F2Physics(frequency_Hz=2400000.));assert np.all(r.query(np.array([0j]))==0)
 assert np.allclose(r.query(np.array([2j]))[0],v[1]*np.exp(1j*np.pi/2*np.arange(1,9)))
 for p in [.5,3.00001]:
  with pytest.raises(ValueError):r.query(np.array([p]))
 with pytest.raises(ValueError):FixedResponse(path,F2Physics(frequency_Hz=2500000.))


def test_archive_outside_observation_not_fabricated():
 from scripts.straight_tube_exp1.f10ac_qd.search import BoundedGridArchive
 a=BoundedGridArchive(solution_dim=2,dims=[2,2],ranges=[[-1,1],[-1,1]])
 result=a.add([[0,0],[1,1]],[5,9],[[0,0],[3,3]])
 assert len(a)==1 and result['status'][1]==0
 assert np.array_equal(a.data()['measures'],[[0,0]])

def test_official_bop_learns_outside_archive_observation():
 from scripts.straight_tube_exp1.f10ac_qd.search import BoundedGridArchive
 from ribs.emitters import BayesianOptimizationEmitter
 from ribs.schedulers import BayesianOptimizationScheduler
 a=BoundedGridArchive(solution_dim=2,dims=[2,2],ranges=[[-1,1],[-1,1]])
 e=BayesianOptimizationEmitter(a,lower_bounds=[0,0],upper_bounds=[1,1],initial_solutions=[[.2,.3],[.6,.7]],seed=45,search_nrestarts=1)
 s=BayesianOptimizationScheduler(a,[e]);s.ask();s.tell([5,9],[[0,0],[3,3]])
 assert len(a)==1 and e.num_evals==2
 assert np.array_equal(e._dataset['measures'],[[0,0],[3,3]])
