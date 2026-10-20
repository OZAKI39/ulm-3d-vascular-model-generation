"""Optuna SQLite ask/evaluate/tell; BoTorch constrained suggestion, native subprocess."""
from common import *
from traction_library import BOUNDS,physical_admission
import optuna,torch,subprocess,time,resource,csv,pickle
from optuna_integration.botorch import BoTorchSampler
torch.set_num_threads(1)
STUDY='stokes_small_native_v1';DB=O/'OPTIMIZATION_STUDY.db'
def candidates(train_x,train_obj,train_con,bounds,pending_x):
 from botorch.models import SingleTaskGP
 from botorch.models.transforms.outcome import Standardize
 from botorch.fit import fit_gpytorch_mll
 from gpytorch.mlls import ExactMarginalLogLikelihood
 from botorch.acquisition.logei import qLogExpectedImprovement
 from botorch.acquisition.objective import GenericMCObjective
 from botorch.sampling.normal import SobolQMCNormalSampler
 from botorch.optim import optimize_acqf
 from botorch.utils.transforms import normalize,unnormalize
 start=time.perf_counter();x=normalize(train_x,bounds);y=torch.cat([train_obj,train_con],dim=-1);gp=SingleTaskGP(x,y,outcome_transform=Standardize(m=y.shape[-1]));fit_gpytorch_mll(ExactMarginalLogLikelihood(gp.likelihood,gp))
 feasible=(train_con<=0).all(-1);assert feasible.any()
 acq=qLogExpectedImprovement(gp,best_f=train_obj[feasible].max(),sampler=SobolQMCNormalSampler(torch.Size([128]),seed=1016),objective=GenericMCObjective(lambda Z,X:Z[...,0]),constraints=[lambda Z,j=j:Z[...,j+1] for j in range(train_con.shape[-1])])
 box=torch.stack([torch.zeros(bounds.shape[-1],dtype=x.dtype),torch.ones(bounds.shape[-1],dtype=x.dtype)])
 candidate,value=optimize_acqf(acq,box,q=1,num_restarts=5,raw_samples=128,options={'batch_limit':5,'maxiter':100});result=unnormalize(candidate.detach(),bounds)
 dump(O/f'BOTORCH_PROPOSAL_after_{len(x):03d}.json',dict(status='ACTUALLY_EXECUTED_EXPLORATORY_ONLY',samples=len(x),train_parameters=train_x.tolist(),train_objectives=train_obj.tolist(),train_constraints=train_con.tolist(),bounds=bounds.tolist(),candidate=result.tolist(),log_expected_improvement=float(value),wall_s=time.perf_counter()-start,scope='Too few points for reliable GP or global optimum; interface/feedback demonstration only',acquisition='qLogExpectedImprovement with two actual outcome constraints'))
 return result
def sampler():
 p=O/'SAMPLER_STATE.pkl'
 if p.exists():
  with p.open('rb') as f:return pickle.load(f)
 return BoTorchSampler(candidates_func=candidates,n_startup_trials=3,independent_sampler=optuna.samplers.RandomSampler(seed=20261016),seed=20261016,device=torch.device('cpu'))
def save_sampler(s):
 p=O/'SAMPLER_STATE.pkl';tmp=p.with_suffix('.tmp')
 with tmp.open('wb') as f:pickle.dump(s,f)
 os.replace(tmp,p)
def export(study):
 rows=[]
 for t in study.trials:
  rows.append(dict(trial=t.number,state=t.state.name,value=t.value,**t.params,**{k:v for k,v in t.user_attrs.items() if not isinstance(v,(list,dict))},**{'constraint_'+k:v for k,v in t.constraints.items()}))
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with (O/'OPTIMIZATION_TRIALS.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def execute_trial(study,trial):
 params={k:trial.suggest_float(k,*v) for k,v in BOUNDS.items()};case=O/'trials'/f'trial_{trial.number:03d}';case.mkdir(parents=True,exist_ok=True)
 completed=sum(t.state==optuna.trial.TrialState.COMPLETE for t in study.trials)
 trial.set_user_attr('phase','INITIAL_RANDOM_DESIGN' if completed<3 else 'BOTORCH_EXPLORATORY_SMALL_DATA');trial.set_user_attr('fidelity','LOW_FIDELITY_PROXY');trial.set_user_attr('case',str(case.relative_to(R)))
 admission=physical_admission(params);dump(case/'PHYSICAL_ADMISSION.json',admission)
 if not admission['valid']:
  trial.set_user_attr('result_class','OUTSIDE_VALIDATED_RANGE');study.tell(trial,state=optuna.trial.TrialState.PRUNED);save_sampler(study.sampler);return
 req=case/'submission.json';dump(req,dict(case=str(case),params=params));save_sampler(study.sampler)
 # Native subprocess keeps its validated Python/MPI environment independent.
 cmd=['/bin/bash','-c','source /workspace/y_clot/env.sh; export LD_LIBRARY_PATH="'+str(LIB.parent)+':$LD_LIBRARY_PATH"; exec /workspace/y_clot/envs/clot/bin/python "'+str(A/'mechanical_evaluator.py')+'" --request "$1"','native',str(req)]
 start=time.perf_counter()
 with (case/'process.log').open('a') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 trial.set_user_attr('subprocess_exit_code',p.returncode);trial.set_user_attr('subprocess_wall_s',time.perf_counter()-start)
 if p.returncode or not (case/'SUMMARY.json').exists():
  trial.set_user_attr('result_class','NUMERICAL_FAILURE');study.tell(trial,state=optuna.trial.TrialState.FAIL);return
 z=json.loads((case/'SUMMARY.json').read_text());row=z['last_row'];trial.set_user_attr('result_class','REAL_NO_FRACTURE' if row['new_softened']==0 else 'REAL_NATIVE_DAMAGE');trial.set_user_attr('history_max_free_nodes',z['history_max_free_nodes']);trial.set_user_attr('W_J',row['W_J']);trial.set_user_attr('D_J',row['D_J']);trial.set_user_attr('wall_s',z['wall_s']);trial.set_user_attr('cpu_s',z['cpu_s']);trial.set_user_attr('peak_RSS_bytes',z['peak_RSS_bytes']);trial.set_user_attr('max_relative_energy_residual',z['max_relative_energy_residual'])
 trial.set_constraint('Re_margin',params['drive_speed_m_s']*25.-.1);trial.set_constraint('energy_residual_margin',z['max_relative_energy_residual']-.01)
 study.tell(trial,z['objective']);save_sampler(study.sampler)
def run(limit=4,propose=False):
 start=time.perf_counter();cpu=time.process_time();s=sampler();study=optuna.create_study(study_name=STUDY,storage='sqlite:///'+str(DB),direction='maximize',sampler=s,load_if_exists=True)
 # Reattach persisted RUNNING trials; never silently assign zero or reset DB.
 pending=[t for t in study.trials if t.state==optuna.trial.TrialState.RUNNING]
 for frozen in pending:execute_trial(study,optuna.trial.Trial(study,frozen._trial_id));export(study)
 while sum(t.state==optuna.trial.TrialState.COMPLETE for t in study.trials)<limit:
  trial=study.ask();execute_trial(study,trial);export(study)
  if study.trials[-1].state==optuna.trial.TrialState.FAIL:break
 restored=optuna.load_study(study_name=STUDY,storage='sqlite:///'+str(DB),sampler=sampler());assert [(t.number,t.state,t.value) for t in restored.trials]==[(t.number,t.state,t.value) for t in study.trials]
 dump(O/f'OPTIMIZER_STAGE_{len(study.trials)}.json',dict(status='PERSISTED_AND_RELOADED',trials=len(study.trials),completed=sum(t.state==optuna.trial.TrialState.COMPLETE for t in study.trials),wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,resume_verified='Fresh Study and sampler deserialized from disk; all states and observations match'))
 if propose and sum(t.state==optuna.trial.TrialState.COMPLETE for t in study.trials)>=limit:
  done=study.get_trials(states=(optuna.trial.TrialState.COMPLETE,));names=sorted(BOUNDS);x=torch.tensor([[t.params[k] for k in names] for t in done],dtype=torch.double);y=torch.tensor([[t.value] for t in done],dtype=torch.double);con=torch.tensor([list(t.constraints.values()) for t in done],dtype=torch.double);bounds=torch.tensor([[BOUNDS[k][j] for k in names] for j in [0,1]],dtype=torch.double)
  torch.manual_seed(20261020);p=candidates(x,y,con,bounds,None)[0];dump(O/'NEXT_CANDIDATE.json',dict(status='PROPOSED_NOT_RUN',params=dict(zip(names,p.tolist())),scope='Small specimen only; no V3 admission without crack evidence'))
 return study
