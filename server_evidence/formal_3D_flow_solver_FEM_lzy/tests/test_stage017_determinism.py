from stage017_helpers import *
import pytest
def test_independent_optimizer_replay_is_identical():
 r=read(OUT/'determinism_comparison.json')
 assert r['status']=='PASS' and all(r['checks'].values())
 assert r['new_volume_meshes_total']<=POLICY['limits']['maximum_total_volume_meshes']
 assert read(OUT/'determinism/optimizer_result.json')['selected_iteration']==result()['selected_iteration']
