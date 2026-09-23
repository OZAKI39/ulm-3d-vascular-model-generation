"""Independent small algebra check; not a substitute for the PDE benchmark."""
import numpy as np


def test_one_global_scalar_constraint_block_and_rhs_sign():
    A=np.diag([2.,3.,4.])
    B=np.array([[1.,-1.,0.]])
    C=np.array([[-1.],[0.],[-1.]])
    matrix=np.block([[A,B.T,C],[B,np.zeros((1,1)),np.zeros((1,1))],
                     [C.T,np.zeros((1,1)),np.zeros((1,1))]])
    assert matrix.shape==(5,5) and C.shape==(3,1) and np.linalg.norm(C)>0
    rhs=np.array([0.,0.,0.,0.,-2.])
    solution=np.linalg.solve(matrix,rhs)
    assert np.linalg.matrix_rank(matrix)==5
    assert abs(float(C[:,0]@solution[:3])+2.)<1e-14
    assert np.linalg.norm(B@solution[:3])<1e-14
    assert solution[-1]>0
    assert np.linalg.norm(matrix@solution-rhs)<1e-14
    np.testing.assert_allclose(np.linalg.solve(matrix,-rhs),-solution,atol=1e-14,rtol=0)
