"""Original frozen V2 error decomposition and classification thresholds."""
import numpy as np

def magnitude(x):return np.linalg.norm(x.reshape(len(x),-1),axis=1)

def stats(x):
 a=magnitude(x);l2=float(np.linalg.norm(a));rms=l2/np.sqrt(len(a))
 return dict(L2=l2,RMS=rms,median=float(np.median(a)),P95=float(np.percentile(a,95)),max=float(a.max()))

def decompose(actual,space,exact):
 total=actual-exact;es=space-exact;ed=actual-space;den=float(np.linalg.norm(exact));ns=float(np.linalg.norm(es));nd=float(np.linalg.norm(ed));nt=float(np.linalg.norm(total))
 cross=float(2*np.sum(es*ed)/den**2)
 point_res=float(np.max(np.abs(total-es-ed)))
 energy_res=abs(nt**2-ns**2-nd**2-2*float(np.sum(es*ed)))
 scale=max(float(np.max(np.abs(actual))),float(np.max(np.abs(space))),float(np.max(np.abs(exact))),np.finfo(float).tiny)
 assert point_res<=32*np.finfo(float).eps*scale
 assert energy_res<=128*np.finfo(float).eps*max(den**2,nt**2,ns**2,nd**2)
 return dict(total_relative_L2=nt/den,space_relative_L2=ns/den,discrete_relative_L2=nd/den,cross_term_normalized=cross,
  alignment_cosine=float(np.sum(es*ed)/(ns*nd)) if min(ns,nd)>64*np.finfo(float).eps*den else None,
  point_identity_max_abs=point_res,energy_identity_abs_residual=energy_res,exact_L2=den,exact_RMS=den/np.sqrt(len(exact)),
  total_absolute=stats(total),space_absolute=stats(es),discrete_absolute=stats(ed))

def classification(d):
 # PASS takes precedence in the overlapping user rules: a passing CFD band remains PASS.
 if d['total_relative_L2']<=.05:return 'PASS'
 if d['space_relative_L2']>.05:
  return 'P1_SPACE_LIMIT_SUPPORTED' if d['discrete_relative_L2']<=.05 else 'MIXED_ERROR'
 return 'DISCRETE_SOLUTION_LIMIT_SUPPORTED'
