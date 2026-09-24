"""Deterministic positive-cone redundancy at the dual Gram matrix resolution.

The dual uses W W.T, so its resolvable singular-value floor is sqrt(n*eps)*||W||.
Linear dependence alone is insufficient for inequalities: a discarded row must
also be represented by NONNEGATIVE coefficients of retained rows. All original
constraints are checked after solving, with the retained system's strict budget.
"""
import numpy as np
from scipy.optimize import nnls


def independent_contact_rows(w,ids):
    singular=np.linalg.svd(w,compute_uv=False)
    tol=np.sqrt(max(w.shape)*np.finfo(float).eps)*singular[0]
    kept=[];dropped=[]
    for k,row in enumerate(w):
        if kept:
            coefficients,_=nnls(w[kept].T,row)
            residual=float(np.linalg.norm(row-coefficients@w[kept]))
            if residual<=tol:
                dropped.append(dict(dropped_constraint_id=ids[k],
                    kept_constraint_ids=[ids[i] for i,c in zip(kept,coefficients) if c>0],
                    nonnegative_coefficients=coefficients.tolist(),residual_norm=residual))
                continue
        kept.append(k)
    after=np.linalg.svd(w[kept],compute_uv=False)
    return kept,dict(method='SVD_GRAM_RESOLUTION_AND_NONNEGATIVE_CONE',
        tolerance=tol,tolerance_formula='sqrt(max(W.shape)*eps)*norm(W,2)',
        singular_values_before=singular.tolist(),singular_values_after=after.tolist(),
        rank_before=int(np.sum(singular>tol)),rank_after=int(np.sum(after>tol)),
        input_constraint_count=len(ids),retained_constraint_count=len(kept),
        kept_constraint_ids=[ids[k] for k in kept],dropped=dropped)
