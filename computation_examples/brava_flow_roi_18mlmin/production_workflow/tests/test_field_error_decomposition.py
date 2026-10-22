"""Vector decomposition must retain cross terms and the original PASS category."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'gmsh_reference'))
from field_error_decomposition import decompose,classification

def test_opposed_errors_cancel_without_subtracting_norms():
 exact=np.array([[0.,0.,1.]])
 p1=exact+np.array([[.1,0.,0.]])
 cfd=exact+np.array([[.02,0.,0.]])
 result=decompose(cfd,p1,exact)
 assert np.isclose(result['cross_term_normalized'],-.016)
 assert np.isclose(result['alignment_cosine'],-1.)
 assert np.isclose(result['total_relative_L2'],.02)
 assert classification(result)=='PASS'

def test_quadratic_nodal_interpolant_is_not_exact_derivative_recovery():
 points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
 nodal=np.zeros((4,3));nodal[:,2]=2-points[:,0]**2-points[:,1]**2
 weights=np.full(4,.25);x=weights@points
 interp=weights@nodal;exact=2-x[0]**2-x[1]**2
 gradient=np.linalg.solve(points[1:]-points[:1],nodal[1:]-nodal[:1]).T
 assert interp[2]==1.5 and exact==1.875
 np.testing.assert_array_equal(gradient[2],[-1.,-1.,0.])
 variance=np.sum(weights*np.sum((points[:,:2]-x[:2])**2,axis=1))
 assert interp[2]-exact==-variance
