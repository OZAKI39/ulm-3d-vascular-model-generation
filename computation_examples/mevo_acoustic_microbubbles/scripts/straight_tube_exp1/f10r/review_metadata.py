"""Repair derived source labels from immutable INPUT, keeping raw records untouched."""
import time
from scripts.straight_tube_exp1.f10p.analyze import layered
from .common import BASE,Q,policy,read,write,identity


def main():
    begin=time.perf_counter();an=BASE/'analysis';runs=read(an/'RUN_ANALYSIS.json')
    gates=read(an/'FINAL_GATE.json');before=gates['groups']['OFF_DT1'].copy();changes={}
    assert not (an/'SOURCE_METADATA_REVIEW.json').exists()
    for path in BASE.glob('*/runs/*/INPUT.json'):
        name=path.parent.name;source=read(path)['source'];original=runs[name]['source']
        runs[name]['original_source_metadata']=original;runs[name]['source']=source
        changes[name]=dict(original_recorded_source=original,input_source=source)
    old=read(Q/'analysis/RUN_ANALYSIS.json');comp=read(an/'COMPARISONS.json')
    off=runs['OFF_DT1'];tops=[off['topology']]+[old['OFF_'+s]['topology'] for s in ['DT05','FINE']]
    update=layered(off,comp['OFF_DT1__DT05'],comp['OFF_DT1__FINE'],tops,policy())
    gates['groups']['OFF_DT1'].update(update)
    write(an/'RUN_ANALYSIS.json',runs);write(an/'FINAL_GATE.json',gates)
    elapsed=time.perf_counter()-begin
    write(an/'SOURCE_METADATA_REVIEW.json',dict(identity=identity(),wall_s=elapsed,changes=changes,
          root_cause='Output field-loop key shadowed schedule-key metadata after integration. All controls, actual forces, input identities and trajectories were correct. Raw RESULT/remote hashes remain unchanged.',
          action='Use immutable INPUT.source for derived analysis; new runner saves source_key before output loops. Restore inherited OFF null-reference rule without altering a threshold or numerical result.',
          before_OFF_gate=before,after_OFF_gate=gates['groups']['OFF_DT1']))
    cost=read(an/'COSTS.json');cost['metadata_review_s']=elapsed;cost['total_scientific_s']+=elapsed
    write(an/'COSTS.json',cost)
    print('Reviewed source metadata in',len(changes),'records; OFF COM:',update['GROUP_POSITION_NUMERICALLY_ADMITTED'])


if __name__=='__main__':main()
