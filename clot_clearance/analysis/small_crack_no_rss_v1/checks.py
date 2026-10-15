from native import *
def topology(ns):
 nb,raw=ns.l.gather_bonds();b=np.ctypeslib.as_array(raw,shape=(3*nb,)).reshape(-1,3)[:,1:].copy();b.sort(axis=1);codes=np.sort(b[:,0]*ns.N+b[:,1])
 fn=ns.lib.lps_event_special;fn.argtypes=[ctypes.c_void_p,ctypes.POINTER(ctypes.c_double)];fn.restype=ctypes.c_int
 count=fn(ns.l.lmp,None);out=np.empty(count);assert fn(ns.l.lmp,ns.ptr(out))==count
 records=[];at=0
 while at<count:
  atom_id=int(out[at]);counts=out[at+1:at+4].astype(int);at+=4;v=out[at:at+counts[2]].astype('<i8');at+=counts[2]
  pieces=[np.sort(v[:counts[0]]),np.sort(v[counts[0]:counts[1]]),np.sort(v[counts[1]:])];digest=hashlib.sha256(np.r_[counts,*pieces].astype('<i8').tobytes()).hexdigest();records.append((atom_id,digest))
 records.sort();return codes,hashlib.sha256(json.dumps(records).encode()).hexdigest()
def compare(p,s):
 err={k:float(np.max(abs(np.asarray(p[k])-np.asarray(s[k])))) for k in ['x','velocity','force','material_force','kappa','d','time','step','native_bonds','Uc','K','native_U']};err['physical_scalars']=float(abs(p['scalars'][:6]-s['scalars'][:6]).max())
 ok=all(err[k]==0 for k in ['x','velocity','kappa','d','time','step','native_bonds']) and max(err['force'],err['material_force'])<=2e-17 and max(err['Uc'],err['K'],err['native_U'],err['physical_scalars'])<=1e-20
 return dict(status='PASS' if ok else 'FAIL',max_absolute_errors=err)
def predict_x(s,h,ref,ids):
 x=s['x'].copy();mask=ref['free_dof'];x[mask]+=h*(s['velocity'][mask]+h/(2*RHO*V)*s['force'][mask]);u,_=boundary(s['time']+h)
 for key,sign in [('bottom_ids',-1),('top_ids',1)]:
  ii=ids[key]-1;x[ii,2]=ref['X'][ii,2]+sign*u
 return x
