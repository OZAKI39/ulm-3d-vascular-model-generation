"""Single deterministic elevation; no mesher, refinement or CFD."""
from common import *
import importlib.util
import numpy as np
import math
spec=importlib.util.spec_from_file_location('frozen_sampler',NATIVE/'p2_sampler.py');sampler=importlib.util.module_from_spec(spec);spec.loader.exec_module(sampler)

def main():
 p=NATIVE/'runtime/PIPE_AFFINE_ELEVATION.npz';assert sha(p)=='a50b1fb8ea5b667947cb4e79b62f16af40008d1bbf83358ae1704d1645877056'
 data=dict(np.load(p));X=data['points_m'].copy();C=data['tetra10'];B=data['boundary6'];tags=data['facet_tags'];edges=data['unique_edges']
 plan=read(BASE/'production_workflow/gmsh_reference/GMSH_REFERENCE_MESH_V2_PLAN.json');R=plan['radius_m'];L=plan['length_m']
 raw=np.load(BASE/'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu/SV_MESH/mesh_arrays.npz');corners=raw['points_m'];n=len(corners)
 wallnodes=np.unique(B[tags==1]);wallmids=wallnodes[wallnodes>=n];r=np.linalg.norm(X[wallmids,:2],axis=1);X[wallmids,:2]*=(R/r)[:,None]
 # J(xi) is affine. Convexity of ||J0^-1 (J-J0)||2 bounds every point
 # in the reference tetrahedron by its four corner values. Bound <1 with
 # detJ0>0 is a whole-reference-domain positivity certificate.
 refs=np.vstack([np.eye(3),np.zeros(3)]);D=sampler.basis(refs)[1]
 bounds=np.empty(len(C));mincorner=np.empty(len(C));maxcorner=np.empty(len(C));det0=np.empty(len(C))
 for start in range(0,len(C),10000):
  sl=slice(start,start+10000);Y=X[C[sl]];Y=Y-Y[:,3,None];J0=Y[:,:3].transpose(0,2,1);inv=np.linalg.inv(J0);det0[sl]=np.linalg.det(J0)
  J=np.einsum('nai,vaj->nvij',Y,D);det=np.linalg.det(J);mincorner[sl]=det.min(1);maxcorner[sl]=det.max(1)
  rel=np.einsum('nij,nvjk->nvik',inv,J-J0[:,None]);bounds[sl]=np.linalg.svd(rel,compute_uv=False)[...,0].max(1)
 certified=(bounds<1)&(det0>0)&(mincorner>0)&np.isfinite(bounds)
 checks=dict(corner_nodes_bitwise_unchanged=np.array_equal(X[:n],corners),topology_unchanged=np.array_equal(np.sort(C[:,:4],axis=1),np.sort(raw['tetra'],axis=1)),edge_midpoint_ID_unique=len(edges)==len(np.unique(edges,axis=0)),shared_midpoint_conformity=True,boundary_labels_unchanged=np.array_equal(tags,raw['facet_tags']),boundary_corner_connectivity_unchanged=np.array_equal(B[:,:3],raw['boundary_triangles']),all_reference_domains_positive_certified=bool(certified.all()),no_duplicate_coordinates=len(np.unique(X,axis=0))==len(X))
 # Every local edge must refer to its globally sorted unique-edge endpoints.
 mapped=edges[C[:,4:]-n];expected=np.sort(C[:,:4][:,sampler.EDGES],axis=-1);checks['shared_midpoint_conformity']=np.array_equal(mapped,expected)
 wallerr=float(np.abs(np.linalg.norm(X[wallnodes,:2],axis=1)-R).max());caperr=max(float(abs(X[np.unique(B[tags==2]),2]).max()),float(abs(X[np.unique(B[tags==3]),2]-L).max()))
 checks['wall_nodes_on_cylinder']=wallerr<=1e-12;checks['cap_planes']=caperr<=plan['gates']['cap_plane_tolerance_m']
 # Exact degree-four Bernstein bound of x^2+y^2-R^2 over every TRI6 wall.
 multi=np.array([(i,j,4-i-j) for i in range(5) for j in range(5-i)])
 lam=multi/4.;bern=np.array([math.factorial(4)/math.prod(math.factorial(int(k)) for k in a)*np.prod(lam**a,axis=1) for a in multi]).T
 shape=np.column_stack([lam*(2*lam-1),4*lam[:,0]*lam[:,1],4*lam[:,1]*lam[:,2],4*lam[:,2]*lam[:,0]])
 world=np.einsum('qa,faj->fqj',shape,X[B[tags==1]])
 nodal_defect=(world[:,:,:2]**2).sum(2)-R**2;coeff=nodal_defect@np.linalg.inv(bern).T
 lo=coeff.min(1);hi=coeff.max(1)
 wall_surface_bound=float(max(abs(np.sqrt(R*R+lo)-R).max(),abs(np.sqrt(R*R+hi)-R).max()))
 checks['whole_wall_radial_Bernstein_bound']=wall_surface_bound<=plan['wall_geometry_limit_m']

 data['points_m']=X;out=HERE/'runtime/PIPE_CURVED_ELEVATION.npz';np.savez_compressed(out,**data)
 np.savez_compressed(HERE/'P2_MESH_QA_SAMPLES.npz',cell_id=np.arange(len(C)),perturbation_bound=bounds,det_affine=det0,det_corner_min=mincorner,det_corner_max=maxcorner,certified=certified)
 result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks={k:bool(v) for k,v in checks.items()},cells=len(C),velocity_nodes=len(X),pressure_corner_nodes=n,unique_edges=len(edges),wall_projected_midpoints=len(wallmids),certified_count=int(certified.sum()),uncertified_cell_ids=np.flatnonzero(~certified).tolist(),negative_corner_cell_ids=np.flatnonzero(mincorner<=0).tolist(),max_perturbation_bound=float(bounds.max()),min_corner_det=float(mincorner.min()),wall_node_error_m=wallerr,wall_surface_radial_error_bound_m=wall_surface_bound,wall_surface_certificate='Degree-four triangular Bernstein convex-hull bounds for x^2+y^2-R^2 on all TRI6 wall facets, FP64 polynomial evaluation',cap_plane_error_m=caperr,certificate='det J0>0 and max at four reference vertices spectral norm(J0^-1 (J-J0))<1; convexity + homotopy certify full element, not quadrature sampling',curved_mesh_sha256=sha(out),affine_input_sha256=sha(p),no_refinement=True,no_corner_movement=True)
 put('P2_CURVED_MESH_QA.json',result);print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
