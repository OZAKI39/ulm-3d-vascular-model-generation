from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_stationary_requires_contact_support(tmp_path):
 from particle_3d.formal_dynamics_p9a5 import AuditSink
 m=dict(completed=False,end_reason='INTEGRATION_SAFETY_STOP',failure_detail='ROUNDOFF_SCALE_STAGNATION')
 assert terminal_status(m,False,12.)=='SOLVER_FAILURE'
 sink=AuditSink(tmp_path/'audit.gz')
 row=dict(velocity=[0.,0.,0.,0.,0.,0.],solver=dict(contact_redundancy={'rank_after':3},contact_kkt={'velocity_budget_m_s':1e-12},multipliers=[1,1,1]))
 for _ in range(31):sink.last.append(deepcopy(row))
 assert not sink.supported()
 sink.last.append(deepcopy(row));assert sink.supported()
 assert terminal_status(m,sink.supported(),12.)=='SUPPORTED_STATIONARY'
 sink.last[-1]['solver']['contact_redundancy']['rank_after']=2;assert not sink.supported()
 sink.close()
