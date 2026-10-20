"""Objective kinematic diagnostics and exact full endpoint-star virtual removal.
No history update, force modification, or new failure threshold is performed.
"""
from common import *
import time

def corotated(ref,x,d):
 b=ref['bonds'];X=ref['X'];n=len(X);nb=len(b);i,j=b.T
 xi=X[j]-X[i];y=x[j]-x[i];w=ref['taper']*(1-d)/ref['r0']**2
 A=np.zeros((n,3,3));C=np.zeros_like(A)
 for a in range(3):
  for c in range(3):
   A[:,a,c]=np.bincount(b.ravel(),weights=np.repeat(w*y[:,a]*xi[:,c],2),minlength=n)
   C[:,a,c]=np.bincount(b.ravel(),weights=np.repeat(w*xi[:,a]*xi[:,c],2),minlength=n)
 ev=np.linalg.eigvalsh(C);valid=ev[:,0]>np.maximum(ev[:,-1],1)*1e-10
 F=A@np.linalg.pinv(C,rcond=1e-10);u,ss,v=np.linalg.svd(F);sgn=np.linalg.det(u@v);u[:,:,2]*=sgn[:,None];rot=u@v
 un=np.einsum('nij,nj->ni',rot[i].transpose(0,2,1),y)-xi;vn=np.einsum('nij,nj->ni',rot[j].transpose(0,2,1),y)-xi;direction=xi/ref['r0'][:,None]
 a=np.sum(un*direction,axis=1);c=np.sum(vn*direction,axis=1)
 shear=.5*(np.linalg.norm(un-a[:,None]*direction,axis=1)+np.linalg.norm(vn-c[:,None]*direction,axis=1))/ref['r0']
 longitudinal=.5*(a+c)/ref['r0'];ok=valid[i]&valid[j]
 return dict(shear=shear,longitudinal=longitudinal,valid=ok,detF=np.linalg.det(F),condition=np.divide(ev[:,-1],np.maximum(ev[:,0],1e-30)))

def virtual_release(model,x,d):
 # U_i = 1/2 sum_a D_a w_a e_a^2 + 1/2 n_i^2/a_i.
 # Failing edge a changes w_a->0 and h_a->1/D_a at both stars;
 # re-minimizing those stars includes every coupled original neighbor.
 e,_=model.geometry(x);ee=e[model.edge];dd=d[model.edge];den=model.Ddir*dd+model.krefdir*(1-dd)
 w=model.krefdir*(1-dd)/den;h=dd/den
 num=np.bincount(model.center,weights=model.l*w*ee,minlength=model.n);a=1/model.A+np.bincount(model.center,weights=model.l**2*h,minlength=model.n)
 nn=num[model.center];aa=a[model.center];newn=nn-model.l*w*ee;newa=aa+model.l**2*(1/model.Ddir-h)
 du=.5*model.Ddir*w*ee**2+.5*(nn**2/aa-newn**2/newa)
 return du[:model.nb]+du[model.nb:]

def main():
 start=time.perf_counter();ref=reference();m=energy_model(ref);bt=ref['bond_types'];summaries=[];checks=[];cands=json.loads((O/'V3_CHUNK_CANDIDATES.json').read_text())['candidates'];rng=np.random.default_rng(20261015)
 for p,s in keyframes():
  step=int(s['step']);g=corotated(ref,s['x'],s['d']);du=virtual_release(m,s['x'],s['d']);base=m.evaluate(s['x'],s['d']);live=(s['d']<1)&(bt<3);below=s['kappa']<S0[bt];sel=live&g['valid']&below&(g['shear']>=.15)
  typed=[]
  for t in [1,2]:
   q=live&(bt==t)&g['valid'];z=q&sel
   typed.append(dict(bond_type=t,live_valid=int(q.sum()),live_invalid=int(np.sum(live&(bt==t)&~g['valid'])),below_history_onset_and_shear_ge_015=int(z.sum()),shear_quantiles=np.quantile(g['shear'][q],[0,.5,.9,.99,1]),virtual_energy_J_quantiles=np.quantile(du[q],[0,.5,.9,.99,1]),max_shear_bond_id=int(np.flatnonzero(q)[np.argmax(g['shear'][q])]+1)))
  examples=np.flatnonzero(sel);examples=examples[np.argsort(du[examples])[-20:][::-1]]
  rec=[]
  for q in examples:rec.append(dict(bond_id=int(q+1),atom_ids=ref['bonds'][q]+1,bond_type=int(bt[q]),stretch=float(base['e'][q]/ref['r0'][q]),corotated_shear=float(g['shear'][q]),longitudinal=float(g['longitudinal'][q]),kappa=float(s['kappa'][q]),d=float(s['d'][q]),virtual_release_J=float(du[q])))
  if step:
   sample=rng.choice(np.flatnonzero(live),8,replace=False)
   err=[]
   for q in sample:
    changed=s['d'].copy();changed[q]=1;exact=base['U']-m.evaluate(s['x'],changed)['U'];err.append(abs(exact-du[q]))
   ferr=np.max(np.abs(base['force']-s['material_force']));checks.append(dict(step=step,full_energy_evaluations=len(sample),virtual_max_abs_error_J=max(err),native_force_max_abs_error_N=ferr,native_energy_abs_error_J=abs(base['U']-float(s['native_U']))))
   assert max(err)<1e-18 and ferr<1e-13
  np.savez_compressed(O/'derived'/f'mixed_step_{step:06d}.npz',bond_ids=np.arange(1,len(bt)+1),**g,virtual_release_J=du,stretch=base['e']/ref['r0'])
  bridges=[]
  for c in cands:
   if c['step']!=step:continue
   bids=np.asarray(c['external_bond_ids'])[np.asarray(c['external_alive'])]-1
   for q in bids:bridges.append(dict(candidate_id=c['candidate_id'],bond_id=int(q+1),valid=bool(g['valid'][q]),shear=float(g['shear'][q]),stretch=float(base['e'][q]/ref['r0'][q]),kappa=float(s['kappa'][q]),d=float(s['d'][q]),virtual_release_J=float(du[q]),below_history_onset=bool(below[q])))
  summaries.append(dict(step=step,time_s=float(s['time']),typed=typed,largest_virtual_release_shear_examples=rec,candidate_bridges=bridges))
  print('MIXED',step,'below_onset_shear015',int(sel.sum()),flush=True)
 # Translation + finite proper rotation invariance and potential-force gradient.
 angle=.73;Q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
 xx=s['x']@Q.T+np.array([.003,-.002,.001]);gg=corotated(ref,xx,s['d']);good=g['valid']&gg['valid'];obj=max(np.max(abs(g['shear'][good]-gg['shear'][good])),np.max(abs(g['longitudinal'][good]-gg['longitudinal'][good])));assert obj<1e-9
 rigid=corotated(ref,ref['X']@Q.T+.001,np.zeros(len(bt)));assert np.max(rigid['shear'][rigid['valid']])<1e-10
 eps=5e-10;grad=[]
 for i in rng.choice(int(ref['clot_count']),8,replace=False):
  for a in range(3):
   xp=s['x'].copy();xm=xp.copy();xp[i,a]+=eps;xm[i,a]-=eps;fd=-(m.evaluate(xp,s['d'])['U']-m.evaluate(xm,s['d'])['U'])/(2*eps);grad.append(abs(fd-base['force'][i,a]))
 assert max(grad)<1e-10
 dump(O/'MIXED_MODE_DIAGNOSTIC.json',dict(status='PASS_DIAGNOSTIC_ONLY',method='Weighted least-squares local F using original live bonds, polar proper rotation at both endpoints; mean corotated transverse relative displacement divided by r0. Rank-deficient local stars marked NOT MEASURED.',threshold_note='0.15 shear is a diagnostic comparison scale only, NOT a critical skew or material failure criterion. kappa must also remain below original typed onset.',virtual_energy='Exact full endpoint-star re-minimization at fixed x after one original bond d->1, checked against independent complete evaluations. Not Gc, not certified fracture dissipation; do not sum single-bond removals as independent energies.',objectivity_max_abs_error=obj,rigid_motion_max_shear=float(np.max(rigid['shear'][rigid['valid']])),force_gradient_max_abs_error_N=max(grad),independent_native_and_virtual_checks=checks,frames=summaries,interpretation='Shear-dominant kinematics below stretch onset can exist, but their presence does not prove a physical failure blind spot. Simple shear stretches oblique bonds already. A calibrated 3D mixed-mode law and non-affine benchmark are required before changing the native criterion.',wall_s=time.perf_counter()-start))
if __name__=='__main__':main()
