from launch import *
assert read(HERE/'P1_RESTART_CORRESPONDENCE.json')['PASS']
assert not list((HERE/'runtime/p1_candidate_fixed_observer/run').glob('algebra_rank_*.jsonl'))
run('p1_candidate_fixed_observer','parity')
