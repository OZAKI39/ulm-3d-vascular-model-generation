"""Tangent equilibrium predictor from an accepted state, using PETSc KSP."""
import numpy as np
from dolfinx import fem
from dolfinx.fem.petsc import assemble_matrix
from petsc4py import PETSc
from .native import vector

def predict(s,accepted_u,prescribed):
 s.u.x.array[:]=accepted_u;s.update();Kfull=assemble_matrix(fem.form(s.Ju));Kfull.assemble();K=assemble_matrix(fem.form(s.Ju),bcs=s.bcs);K.assemble()
 known=K.createVecRight();known.set(0);known.array[s.fixed_flat]=prescribed[s.fixed_flat]-accepted_u[s.fixed_flat];b=K.createVecLeft();Kfull.mult(known,b);b.scale(-1);b.array[:]-=vector(s.Ru);b.array[s.fixed_flat]=known.array_r[s.fixed_flat]
 du=K.createVecRight();ksp=PETSc.KSP().create(s.domain.comm);ksp.setOperators(K);ksp.setType('preonly');pc=ksp.getPC();pc.setType('lu');pc.setFactorSolverType('mumps');ksp.solve(b,du)
 if ksp.getConvergedReason()<=0:raise RuntimeError('Tangent predictor linear solve failed')
 s.u.x.array[:]=accepted_u+du.array_r;s.u.x.array[s.fixed_flat]=prescribed[s.fixed_flat];s.update();info=dict(J_min=s.minJ,J_max=s.maxJ,du_max_m=float(abs(du.array_r).max()),KSP_reason=int(ksp.getConvergedReason()))
 ksp.destroy();K.destroy();Kfull.destroy();known.destroy();b.destroy();du.destroy();return info
