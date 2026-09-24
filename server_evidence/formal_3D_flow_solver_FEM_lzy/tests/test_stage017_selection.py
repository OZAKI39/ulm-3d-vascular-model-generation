from stage017_helpers import *
import pytest
from fem3d.adaptive_port import select_iteration,choose_surface

def test_select_cheapest_acceptable_history_even_if_quality_less_pretty():
 def row(name,dofs,tet,low):return {'iteration':name,'surface_status':'PASS','acceptance':{'status':'PASS'},'volume':{'proxy':{'N_P2_velocity_proxy':dofs,'N_tetra':tet},'quality':{'cap_adjacent_below_0_1':low}}}
 assert select_iteration([row('dense',100,100,0),row('cheap',90,110,1)])=='cheap'
 assert select_iteration([row('more_cells',90,110,0),row('fewer_cells',90,100,1)])=='fewer_cells'

def test_actual_surface_selection_uses_minimum_tested_triangles():
 for p,s in result()['initial_surface_searches'].items():
  chosen=choose_surface(s['trials']);assert chosen['trial']==s['selected_trial']
 assert select_iteration(result()['iterations'])==result()['selected_iteration']
