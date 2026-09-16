from pathlib import Path
import numpy as np,json,hashlib,subprocess,os,shutil
R=Path(__file__).resolve().parents[1];W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work']);D=W/'scalar_work_beta0125';U=R/'input/pecnut_scalars';beta=.125
near=np.load(D/'scalars_general_resistance_nearfield.npy');mid=np.load(D/'scalars_general_resistance_midfield.npy');new=np.concatenate([near,mid],axis=2);assert new.shape==(11,2,38,2) and np.isfinite(new).all()
oldbase=np.loadtxt(U/'values_of_lambda.txt');oldlambda=sorted(set(oldbase)|set(1/oldbase));newbase=sorted(set(oldbase)|{beta});lambdas=sorted(set(oldlambda)|{beta,1/beta});old=np.load(U/'scalars_general_resistance.npy');merged=np.stack([new[:,:,:,[beta,1/beta].index(v)] if v in [beta,1/beta] else old[:,:,:,oldlambda.index(v)] for v in lambdas],axis=3)
np.save(D/'scalars_general_resistance.npy',merged);(D/'values_of_lambda.txt').write_text(''.join(f'{v:.17g}\n' for v in newbase));(D/'find_resistance_scalars').mkdir(exist_ok=True)
with (R/'logs/SCALAR_SUBTRACT_BETA0125.log').open('w') as f:p=subprocess.run([str(W/'env/bin/python'),'-B',str(D/'subtract_R2Binfinity_from_scalars.py')],cwd=D,env=dict(os.environ,PYTHONPATH=str(W/'runtime/stokesian_dynamics'),OPENBLAS_NUM_THREADS='1'),stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0
delta=np.load(D/'find_resistance_scalars/scalars_general_resistance_d.npy');assert delta.shape==merged.shape and np.isfinite(delta).all()
old_delta=np.load(U/'scalars_general_resistance_d.npy');inherit_error=max(float(np.max(abs(delta[:,:,:,lambdas.index(v)]-old_delta[:,:,:,i]))) for i,v in enumerate(oldlambda));assert inherit_error<=1e-12
for n in ['values_of_lambda.txt','values_of_s_dash.txt','scalars_general_resistance.npy']:shutil.copy2(D/n,U/n)
shutil.copy2(D/'find_resistance_scalars/scalars_general_resistance_d.npy',U/'scalars_general_resistance_d.npy')
for p in U.glob('*'):shutil.copy2(p,W/'runtime/stokesian_dynamics/resistance_scalars'/p.name)
np.savez_compressed(R/'raw/SCALAR_GENERATION_BETA0P125_FIRST.npz',near=near,mid=mid,merged=merged,excess=delta,lambdas=lambdas)
report=dict(status='PASS',generated_lambda=[beta,1/beta],inherited_lambda=oldlambda,inherit_max_absolute_delta_change=inherit_error,grid_points=38,source_math_changes=False,input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in U.glob('*')});(R/'validation/SCALAR_GENERATION_BETA0P125.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
