import datetime,h5py,numpy as np
from scripts.straight_tube_exp1.f10ab.prepare import shape
from scripts.straight_tube_exp1.f10ab.metrics import geometry
from scripts.straight_tube_exp1.f10s.common import context
from .common import *
def main():
 BASE.mkdir(parents=True,exist_ok=True);AD.mkdir(parents=True,exist_ok=True)
 assert not (BASE/'PREREGISTRATION.json').exists()
 cfg=dict(total_scientific_wall_budget_s=5400,per_trajectory_s=120,physical_stage_budget_s=2400,optimizer_comparison_s=1200,extension_s=600,fields_s=300,apecss_s=900,postprocessing_s=1800,
 comparison_evaluations_per_method=32,common_initial=8,winner_additional=32,max_workers=4,parallel_pilot_workers=[1,2,4],comparison_duration_s=.002,dt_s=1e-6,frequencies_Hz=[2400000.,2500000.,2600000.],
 descriptor=['delta_log_Rg_3d','delta_log_aspect_3d'],ranges=[[-.08,.08],[-.16,.16]],grid=[8,8],out_of_range='retained but not archived, never clipped',
 quality='100/(1+||COM drift||/1um); valid numerical full window only; Task A permits splitting; Task B requires endpoint connected, sampled path reported separately',
 selection='highest valid archive coverage per total wall, then QD score per wall; one stochastic run only; no universal winner claim',
 action_space='G2: [common X amplitude 200..1100 Pa-eq, X+ relative phase 150..240 deg, common Y amplitude 0..1100 Pa-eq, Y+ phase 0..360 deg]; X- and Y- phases zero',
 numerical_floor='similarity residual > max(1nm,5*dt/cache uncertainty); descriptor separation > max(1e-4,5*sum uncertainties); three unseen seeds, shared sign required',
 heldout_seeds=[301801,301802,301803],field_pressure_atol=1e-4,field_gradient_rtol=.002,response_harmonic_rtol=.01,response_energy_rtol=.02,flow='STRICT_ZERO',radius_m=R,labels=LABELS,training_ready=False)
 write(CONFIG,cfg)
 with h5py.File(AB/'INITIAL_STATES.h5') as old,h5py.File(BASE/'INITIAL_STATES.h5','x') as h:
  for k in ['S0_r0','S1_r0','S2_r0','S4_r0','S6_r0']:h[k]=old[k][:]
  for i,seed in enumerate(cfg['heldout_seeds']):h['S2_unseen_'+str(i)]=shape('S2',seed)
  audit={}
  c=context()
  for k in h:
   x=h[k][:];audit[k]=geometry(x);assert min(c['cavity'].closest_wall(q)[0]-R for q in x)>100e-6
 write(BASE/'INITIAL_STATE_AUDIT.json',audit)
 protected={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'src/mevo_mb').rglob('*.py')}
 for p in [CONFIG,BASE/'INITIAL_STATES.h5',OLD/'inputs/response.h5',OLD/'cache/coarse.h5',OLD/'cache/fine.h5',AB/'INPUT_OUTPUT_SHA256.json']:
  protected[str(p.relative_to(ROOT))]=sha(p)
 write(BASE/'PREREGISTRATION.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),config=cfg,protected_sha256=protected,identity=identity(),baseline_commit='9d304c22817599f85079ada5a42b2e143f3fe055'))
if __name__=='__main__':main()
