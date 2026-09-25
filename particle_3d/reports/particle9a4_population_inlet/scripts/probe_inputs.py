from pathlib import Path
import importlib.util,time,json,hashlib,sys,numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
R=ROOT/'particle_3d/reports/particle9a4_population_inlet';old=ROOT/'particle_3d/reports/network_derived_flow_mb_validation_v1'
s=importlib.util.spec_from_file_location('paired_runner',old/'scripts/runner.py');runner=importlib.util.module_from_spec(s);s.loader.exec_module(runner)
from particle_3d.inlet_flux import frozen_boundary_flux
from particle_3d.injection_method_c import TruncatedSonoVue
from particle_3d.injection_admission import FiniteSizeAdmission
import pyvista as pv
t=time.perf_counter();bundle=old/'server_bundle';env=runner.make_environment(bundle,'NEW');flow=pv.read(bundle/'inputs/NEW.vtu');flux,samplers=frozen_boundary_flux(env.mesh,flow,env.boundaries);dist=TruncatedSonoVue(ROOT/'sonovue_size_distribution_v0');sampler=samplers['INLET'];checker=FiniteSizeAdmission(None,None,wall=env.wall,field=env.field)
start=time.perf_counter();accepted=0;reasons={}
for i in range(100):
 d,_=dist.sample(np.random.default_rng([2026092594,i,941]));xyz,ids=sampler.sample(np.random.default_rng([2026092594,i,942]));event=dict(species='MB',particle_id=i+1,radius_m=d/2,q=[1,0,0,0]);p,status,detail=checker.check(event,xyz[0],{});accepted+=p is not None;reasons[status]=reasons.get(status,0)+1
result=dict(flow_sha256=runner.sha(bundle/'inputs/NEW.vtu'),flux=flux,F_original_4um=dist.mass,source=dist.contract,probe_accepted=accepted,probe_reasons=reasons,probe_100_runtime_s=time.perf_counter()-start,runtime_s=time.perf_counter()-t,wall_provenance=env.wall.provenance,wall_roundoff_m=env.wall.roundoff_m)
(R/'data/input_probe.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source','wall_provenance']},indent=2))
