import numpy as np
from scipy import sparse
from particle_3d.resistance_assembly import ResistanceSystem
from particle_3d.resistance_solver import solve_resistance
from particle_3d.kinematic_contact import ContactConstraint

def test_actual_id7_duplicate(id7):
 f=id7;sys=ResistanceSystem((7,),sparse.csr_matrix(f['R']),np.array(f['b']),np.array(f['self_diagonal']),np.array(f['free']),[],[])
 cs=[ContactConstraint(tuple(r['canonical_id']),7,None,r['normal'],r['particle_point_m']) for r in f['contacts']]
 solved=solve_resistance(sys,constraints=cs);a=solved.record
 assert np.isfinite(solved.velocity).all()
 assert a['contact_redundancy']['retained_constraint_count']==1
 assert a['contact_redundancy']['rank_before']==a['contact_redundancy']['rank_after']==1
 assert a['contact_kkt']['condition']<2
 assert min(a['normal_speeds_m_s'])>=-a['contact_kkt']['velocity_budget_m_s']
 assert max(abs(np.array(f['J'])@solved.velocity))<1e-15
