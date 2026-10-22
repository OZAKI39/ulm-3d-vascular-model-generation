"""Official optimizer integration on frozen measured trajectories only.
This offline replay is NOT formal physical MF optimization, NOT new trajectories,
and does NOT establish resource savings or screening false-negative probability.
"""
import time,traceback
import numpy as np
from scipy.stats import qmc
from .common import *

def main():
    begin=time.perf_counter();actions=['OFF','X_PAIR','Y_PAIR','XY_INPHASE','XY_QUADRATURE','ROT_PLUS','ROT_MINUS','X_SINGLE'];data={}
    for i,a in enumerate(actions):
        data[i]={}
        for level,sub in [(1,'lf'),(2,'hf')]:
            r=read(BASE/sub/('pair_r0_'+a)/'SUMMARY.json');assert r['status']=='PASS';data[i][level]=dict(value=r['final']['stripe_order']-r['initial']['stripe_order'],measured_wall_s=r.get('task_wall_s',r.get('wall_s')),source=str(BASE/sub/('pair_r0_'+a)/'SUMMARY.json'),task=r['task'])
    scope='OFFLINE_API_REPLAY_NOT_SCIENTIFIC_SEARCH';write(BASE/'OPTIMIZER_REPLAY_PLAN.json',dict(scope=scope,actions=actions,metric='Measured stripe-order gain at 2 ms; not pattern confirmation',same_physical_duration=True,formal_MF_admission=False,reason='No HF-confirmed pattern; false-negative sensitivity unidentifiable',identity=identity()));outputs={}
    # SMAC3's own Hyperband intensifier. Budget denotes integrator fidelity, NOT time.
    try:
        from ConfigSpace import ConfigurationSpace,Categorical
        from smac import MultiFidelityFacade,Scenario
        start=time.perf_counter();trace=[];cs=ConfigurationSpace(seed=462070);cs.add(Categorical('action',list(range(8))))
        def objective(config,seed=0,budget=1):
            i=int(config['action']);level=1 if budget<2 else 2;r=data[i][level];trace.append(dict(action=actions[i],level=level,**{k:v for k,v in r.items() if k!='task'}));return -r['value']
        scenario=Scenario(cs,deterministic=True,n_trials=16,min_budget=1,max_budget=2,seed=462070,walltime_limit=120,output_directory=BASE/'optimizer_replay/smac',name='f10af_frozen_replay')
        smac=MultiFidelityFacade(scenario,objective,intensifier=MultiFidelityFacade.get_intensifier(scenario,eta=2,incumbent_selection='highest_budget'),initial_design=MultiFidelityFacade.get_initial_design(scenario,n_configs=4),overwrite=False)
        inc=smac.optimize();outputs['SMAC3']=dict(status='PASS',actual_overhead_s=time.perf_counter()-start,trace=trace,recommended_action=actions[int(inc['action'])],scope=scope)
    except Exception as ex:outputs['SMAC3']=dict(status='FAIL',reason=repr(ex),traceback=traceback.format_exc(),scope=scope)
    # Official multi-fidelity GP + cost-aware qMFKG; evaluate only measured actions.
    try:
        import torch
        from botorch import fit_gpytorch_mll
        from botorch.models.gp_regression_fidelity import SingleTaskMultiFidelityGP
        from botorch.models.transforms.outcome import Standardize
        from botorch.models.cost import AffineFidelityCostModel
        from botorch.acquisition.cost_aware import InverseCostWeightedUtility
        from botorch.acquisition.knowledge_gradient import qMultiFidelityKnowledgeGradient
        from botorch.acquisition.utils import project_to_target_fidelity
        from gpytorch.mlls import ExactMarginalLogLikelihood
        start=time.perf_counter();torch.set_num_threads(1);torch.manual_seed(462071);dtype=torch.double;features=[]
        for i in range(8):
            t=data[i][2]['task'];features.append(np.r_[np.array(t['amplitudes'])/1500,np.mod(t['phases'],2*np.pi)/(2*np.pi)])
        features=np.asarray(features);seen=[(0,1),(1,1),(3,1),(6,1),(0,2),(1,2)];trace=[];bounds=torch.tensor([[0.]*9,[1.]*9],dtype=dtype)
        lf=float(np.median([data[i][1]['measured_wall_s'] for i in data]));hf=float(np.median([data[i][2]['measured_wall_s'] for i in data]));cost=AffineFidelityCostModel(fidelity_weights={8:max(hf-lf,.001)},fixed_cost=lf)
        for iteration in range(2):
            X=torch.tensor(np.array([np.r_[features[i],l-1] for i,l in seen]),dtype=dtype);Y=torch.tensor([[data[i][l]['value']] for i,l in seen],dtype=dtype);model=SingleTaskMultiFidelityGP(X,Y,train_Yvar=torch.full_like(Y,1e-10),data_fidelities=[8],outcome_transform=Standardize(m=1));fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood,model),optimizer_kwargs={'options':{'maxiter':60}})
            targets=torch.tensor(np.c_[features,np.ones(8)],dtype=dtype)
            with torch.no_grad():current=model.posterior(targets).mean.max()
            acq=qMultiFidelityKnowledgeGradient(model=model,num_fantasies=4,current_value=current,cost_aware_utility=InverseCostWeightedUtility(cost),project=lambda X:project_to_target_fidelity(X,{8:1.},d=9))
            choices=[(i,l) for i in range(8) for l in [1,2] if (i,l) not in seen];values=[]
            for i,l in choices:
                xx=torch.tensor(np.r_[features[i],l-1],dtype=dtype).reshape(1,1,9)
                val=acq.evaluate(xx,bounds=bounds,num_restarts=2,raw_samples=8,options={'maxiter':15});values.append(float(val.detach().reshape(-1)[0]))
            j=int(np.argmax(values));i,l=choices[j];seen.append((i,l));trace.append(dict(action=actions[i],level=l,acquisition=values[j],measured_objective=data[i][l]['value']))
        outputs['BoTorch_qMFKG']=dict(status='PASS',actual_overhead_s=time.perf_counter()-start,trace=trace,initial_measured_observations=6,scope=scope,limitation='Offline acquisition exercise; relaxed inner hypothetical search, no newly executed controls or proof of scientific screening')
    except Exception as ex:outputs['BoTorch_qMFKG']=dict(status='FAIL',reason=repr(ex),traceback=traceback.format_exc(),scope=scope)
    start=time.perf_counter();ids=np.floor(qmc.Sobol(1,scramble=True,seed=462072).random_base2(3).ravel()*8).astype(int);outputs['Sobol_replay']=dict(status='PASS',actual_overhead_s=time.perf_counter()-start,trace=[dict(action=actions[i],level=2,value=data[int(i)][2]['value']) for i in ids],scope=scope)
    write(BASE/'OPTIMIZER_API_REPLAY.json',dict(methods=outputs,formal_MF_optimization='NOT_RUN',new_physics_calls=0,resource_saving_comparison='INCONCLUSIVE; replay overhead cannot establish saved physics compute',CMA_MAE_comparison='NOT_RUN; historical objectives and duration differ',wall_s=time.perf_counter()-begin,identity=identity()))
    print('TOOLING', {k:v['status'] for k,v in outputs.items()},flush=True)
if __name__=='__main__':main()
