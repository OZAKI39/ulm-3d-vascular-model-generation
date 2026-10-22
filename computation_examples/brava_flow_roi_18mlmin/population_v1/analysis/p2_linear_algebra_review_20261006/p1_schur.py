from launch import *
assert read(HERE/'P1_RESTART_CORRESPONDENCE.json')['PASS']
run('p1_schur','parity')
