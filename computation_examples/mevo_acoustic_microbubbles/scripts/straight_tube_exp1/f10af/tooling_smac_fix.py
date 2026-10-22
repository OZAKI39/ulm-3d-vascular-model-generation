"""Bound finite catalog replay to its actual configuration count; preserve v1."""
import time
from ConfigSpace import ConfigurationSpace,Categorical
from smac import MultiFidelityFacade,Scenario
from .common import *
def main():
    old=read(BASE/'OPTIMIZER_API_REPLAY.json');start=time.perf_counter();actions=['OFF','X_PAIR','Y_PAIR','XY_INPHASE','XY_QUADRATURE','ROT_PLUS','ROT_MINUS','X_SINGLE'];trace=[];cs=ConfigurationSpace(seed=462070);cs.add(Categorical('action',list(range(8))))
    def objective(config,seed=0,budget=1):
        i=int(config['action']);level=1 if budget<2 else 2;r=read(BASE/('lf' if level==1 else 'hf')/('pair_r0_'+actions[i])/'SUMMARY.json');value=r['final']['stripe_order']-r['initial']['stripe_order'];trace.append(dict(action=actions[i],level=level,value=value,measured_source_wall_s=r.get('task_wall_s',r['wall_s'])));return -value
    scenario=Scenario(cs,deterministic=True,n_trials=8,min_budget=1,max_budget=2,seed=462070,walltime_limit=120,output_directory=BASE/'optimizer_replay/smac_v2',name='bounded_frozen_replay')
    smac=MultiFidelityFacade(scenario,objective,intensifier=MultiFidelityFacade.get_intensifier(scenario,eta=2,incumbent_selection='highest_budget'),initial_design=MultiFidelityFacade.get_initial_design(scenario,n_configs=4),overwrite=False);inc=smac.optimize();old['v1_sha256']=sha(BASE/'OPTIMIZER_API_REPLAY.json');old['methods']['SMAC3']=dict(status='PASS',actual_overhead_s=time.perf_counter()-start,trace=trace,recommended_action=actions[int(inc['action'])],scope='OFFLINE_API_REPLAY_NOT_SCIENTIFIC_SEARCH',v1_failure='Finite configuration space exhausted with requested 16 trials; preserved unchanged');old['identity']=identity();write(BASE/'OPTIMIZER_API_REPLAY_V2.json',old)
if __name__=='__main__':main()
