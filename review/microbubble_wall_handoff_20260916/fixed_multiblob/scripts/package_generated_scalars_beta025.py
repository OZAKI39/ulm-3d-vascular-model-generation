from pathlib import Path
import numpy as np,json,hashlib,subprocess,os,shutil
R=Path(__file__).resolve().parents[1];W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work']);D=W/'scalar_work_beta025';U=R/'input/pecnut_scalars'
near=np.load(D/'scalars_general_resistance_nearfield.npy');mid=np.load(D/'scalars_general_resistance_midfield.npy');new=np.concatenate([near,mid],axis=2);assert new.shape==(11,2,38,2) and np.isfinite(new).all()
base=np.load(U/'scalars_general_resistance.npy');original_lambdas=np.array([.01,.1,.5,1,2,10,100]);lambdas=np.array([.01,.1,.25,.5,1,2,4,10,100]);merged=np.empty((11,2,38,len(lambdas)))
for i,lam in enumerate(lambdas):merged[:,:,:,i]=new[:,:,:,[.25,4].index(lam)] if lam in [.25,4] else base[:,:,:,list(original_lambdas).index(lam)]
np.save(D/'scalars_general_resistance.npy',merged);(D/'values_of_lambda.txt').write_text('.01\n.1\n.25\n.5\n1\n');(D/'find_resistance_scalars').mkdir(exist_ok=True)
env=dict(os.environ,PYTHONPATH=str(W/'runtime/stokesian_dynamics'),OPENBLAS_NUM_THREADS='1')
with (R/'logs/SCALAR_SUBTRACT_BETA025.log').open('w') as f:p=subprocess.run([str(W/'env/bin/python'),'-B',str(D/'subtract_R2Binfinity_from_scalars.py')],cwd=D,env=env,stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0,(R/'logs/SCALAR_SUBTRACT_BETA025.log').read_text()[-1500:]
delta=np.load(D/'find_resistance_scalars/scalars_general_resistance_d.npy');assert delta.shape==merged.shape
out=R/'input/pecnut_scalars';out.mkdir(exist_ok=True)
for name in ['values_of_lambda.txt','values_of_s_dash.txt','scalars_general_resistance.npy']:shutil.copy2(D/name,out/name)
shutil.copy2(D/'find_resistance_scalars/scalars_general_resistance_d.npy',out/'scalars_general_resistance_d.npy')
for p in out.glob('*'):shutil.copy2(p,W/'runtime/stokesian_dynamics/resistance_scalars'/p.name)
np.savez_compressed(R/'raw/SCALAR_GENERATION_BETA0P25_FIRST.npz',near=near,mid=mid,merged=merged,excess=delta,lambdas=lambdas)
report={'status':'PASS','generated_lambda':[.25,4.],'inherited_lambda':original_lambdas.tolist(),'grid_points':38,'source_math_changes':False,'adapter':'concatenate upstream near/mid arrays; merge unchanged inherited lambda columns; unmodified upstream subtract-R2Binfinity script','upstream_IO_workaround':'combine_scalars.py human table has list.shape error; use direct array concatenation; create expected output subdirectory for subtraction','input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*')}};(R/'validation/SCALAR_GENERATION_BETA0P25.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
