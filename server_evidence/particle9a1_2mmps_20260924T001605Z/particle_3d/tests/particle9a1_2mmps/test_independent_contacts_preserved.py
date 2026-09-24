import numpy as np
from particle_3d.contact_redundancy import independent_contact_rows

def test_distinct_and_opposite_nonpenetration_information():
 for rows in [np.eye(3),np.array([[1,0],[0,1],[-1,-1]]),np.array([[1,0],[-1,0]]),np.array([[1,0],[1,1e-5]])]:
  kept,_=independent_contact_rows(rows,list(range(len(rows))))
  assert kept==list(range(len(rows)))
 kept,a=independent_contact_rows(np.array([[1.,0],[0,1],[1,1]]),[0,1,2])
 assert kept==[0,1] and a['dropped'][0]['dropped_constraint_id']==2
