PetscErrorCode VecGetLocalVector(Vec v, Vec w)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(v, VEC_CLASSID, 1);
  PetscValidHeaderSpecific(w, VEC_CLASSID, 2);
  VecCheckSameLocalSize(v, 1, w, 2);
  if (v->ops->getlocalvector) {
    PetscUseTypeMethod(v, getlocalvector, w);
  } else {
    PetscScalar *a;

    PetscCall(VecGetArray(v, &a));
    PetscCall(VecPlaceArray(w, a));
  }
  PetscCall(PetscObjectStateIncrease((PetscObject)w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecRestoreLocalVector(Vec v, Vec w)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(v, VEC_CLASSID, 1);
  PetscValidHeaderSpecific(w, VEC_CLASSID, 2);
  if (v->ops->restorelocalvector) {
    PetscUseTypeMethod(v, restorelocalvector, w);
  } else {
    PetscScalar *a;
    PetscCall(VecGetArray(w, &a));
    PetscCall(VecRestoreArray(v, &a));
    PetscCall(VecResetArray(w));
  }
  PetscCall(PetscObjectStateIncrease((PetscObject)w));
  PetscCall(PetscObjectStateIncrease((PetscObject)v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecGetArray(Vec x, PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(x, VEC_CLASSID, 1);
  PetscCall(VecSetErrorIfLocked(x, 1));
  if (x->ops->getarray) { /* The if-else order matters! VECNEST, VECCUDA etc should have ops->getarray while VECCUDA etc are petscnative */
    PetscUseTypeMethod(x, getarray, a);
  } else if (x->petscnative) { /* VECSTANDARD */
    *a = *((PetscScalar **)x->data);
  } else SETERRQ(PetscObjectComm((PetscObject)x), PETSC_ERR_SUP, "Cannot get array for vector type \"%s\"", ((PetscObject)x)->type_name);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecRestoreArray(Vec x, PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(x, VEC_CLASSID, 1);
  if (a) PetscValidPointer(a, 2);
  if (x->ops->restorearray) {
    PetscUseTypeMethod(x, restorearray, a);
  } else PetscCheck(x->petscnative, PetscObjectComm((PetscObject)x), PETSC_ERR_SUP, "Cannot restore array for vector type \"%s\"", ((PetscObject)x)->type_name);
  if (a) *a = NULL;
  PetscCall(PetscObjectStateIncrease((PetscObject)x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecGetArrayWrite(Vec x, PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(x, VEC_CLASSID, 1);
  PetscValidPointer(a, 2);
  PetscCall(VecSetErrorIfLocked(x, 1));
  if (x->ops->getarraywrite) {
    PetscUseTypeMethod(x, getarraywrite, a);
  } else {
    PetscCall(VecGetArray(x, a));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecRestoreArrayWrite(Vec x, PetscScalar **a)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(x, VEC_CLASSID, 1);
  if (a) PetscValidPointer(a, 2);
  if (x->ops->restorearraywrite) {
    PetscUseTypeMethod(x, restorearraywrite, a);
  } else if (x->ops->restorearray) {
    PetscUseTypeMethod(x, restorearray, a);
  }
  if (a) *a = NULL;
  PetscCall(PetscObjectStateIncrease((PetscObject)x));
  PetscFunctionReturn(PETSC_SUCCESS);
}