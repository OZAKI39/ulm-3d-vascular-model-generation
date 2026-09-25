import hashlib,json,numpy as np
from .conftest import ROOT
from particle_3d.inlet_flux import InletFluxSampler,positive_pieces
from particle_3d.continuous_infusion import counter_rng,P9A4_POSITION

def test_original_source_hash_and_exact_position(real_source):
 audit=json.loads((ROOT/'particle_3d/reports/particle9a4_population_inlet/data/code_logic_audit.json').read_text())
 r=next(r for r in audit['findings'] if r['file'].endswith('/inlet_flux.py'))
 assert hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest()==r['source_sha256']
 for sid in range(1,21):
  row=real_source.proposal(sid);xyz,tri=real_source.sampler.sample(counter_rng(real_source.master_seed,sid,P9A4_POSITION))
  assert np.array_equal(row['position_m'],xyz[0]) and row['inlet_triangle_id']==tri[0]

def test_zero_backflow_clipping_exact_flux():
 xyz=np.array([[0,0,0],[1,0,0],[0,1,0.]])
 pieces=positive_pieces(xyz,[-1,1,0])
 assert abs(sum(p[2] for p in pieces)-1/12)<1e-15
 sampler=InletFluxSampler(xyz[None],np.array([[0.,1.,2.]]))
 x,_=sampler.sample(np.random.default_rng(8194),40000)
 assert np.max(abs(x.mean(0)-sampler.expectation()))<.006
