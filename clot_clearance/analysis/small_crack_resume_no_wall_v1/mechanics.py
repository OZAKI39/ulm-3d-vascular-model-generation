from native import *
from spatial_model import SpatialLPS
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
class Mechanics:
 def __init__(self,ns,ids,initial,gates):
  self.ns=ns;self.ids=ids;self.ref=ns.ref;self.b=ns.ref['bonds'];self.N=ns.N;self.NB=ns.NB;self.g=gates;self.mask=self.ref['free_dof'];self.bc=np.r_[ids['bottom_ids'],ids['top_ids']]-1;self.initial=initial;self.model=SpatialLPS(self.ref);self.contact=Contact(self.ref);self.W=0.;self.Jc=0.;self.cwork=0.;self.scale=1e-15;self.defects=np.zeros(4);self.Dnode=np.zeros(self.N);self.maxerr={};self.last=self.model.evaluate(initial['x'],initial['d']);self.prevfc=self.contact.evaluate(initial['x'],initial['d'])[1];self.components=1;self.labels=np.zeros(self.N,int);self.update_components(initial)
  self.Xaug=np.c_[self.ref['X'],np.ones(self.N)];self.fit=np.linalg.pinv(self.Xaug)
 def error(self,name,val):self.maxerr[name]=max(self.maxerr.get(name,0.),float(val))
 def update_components(self,s):
  alive=self.b[s['d']<1];graph=coo_matrix((np.ones(len(alive)*2),(alive.T.ravel(),alive[:,::-1].T.ravel())),shape=(self.N,self.N)).tocsr();self.components,self.labels=connected_components(graph,directed=False)
 def check(self,p,s,h):
  g=self.g;b=self.b;mask=self.mask;bc=self.bc;dx=s['x']-p['x'];dv=s['velocity']-p['velocity'];sc=s['scalars'];dd=s['d']-p['d'];newfail=np.flatnonzero((s['d']==1)&(p['d']<1));newsoft=np.flatnonzero((s['d']>0)&(p['d']==0));natural=(s['d']==1)&~self.ref['initial_cut'];new=self.model.evaluate(s['x'],s['d']);old=self.model.evaluate(s['x'],p['d']);kk,dr=self.model.update(s['x'],p['kappa'],s_complete=SC)
  self.error('material_force_N',abs(new['force']-s['material_force']).max());self.error('U_J',abs(new['U']-sc[0]));self.error('D_increment_J',abs(old['U']-new['U']-sc[2]));self.error('kappa',abs(kk-s['kappa']).max());self.error('d',abs(dr-s['d']).max());self.error('history_d_rollback',np.max(-dd));self.error('history_kappa_rollback',np.max(p['kappa']-s['kappa']));self.error('KE_J',abs(.5*RHO*V*np.sum(s['velocity']**2)-s['K']));self.error('failed_tension_N',s['failed_tension_max_N'])
  dnode=old['node_energy']-new['node_energy'];self.Dnode+=dnode;self.error('star_Dsum_J',abs(self.Dnode.sum()-sc[1]));self.error('star_U_sum_J',abs(new['node_energy'].sum()-sc[0]));self.error('negative_star_D_J',np.max(-dnode))
  pairs=self.contact.tree(s['x']).query_pairs(DX,output_type='ndarray');pairs=pairs[~np.isin(pairs[:,0]*self.N+pairs[:,1],self.contact.codes[s['d']<1])];penetration=0.;relative_contact_step=0.
  if len(pairs):
   dist=np.linalg.norm(s['x'][pairs[:,1]]-s['x'][pairs[:,0]],axis=1);penetration=float(np.max(1-dist/DX));relative_contact_step=float(np.linalg.norm(dx[pairs[:,1]]-dx[pairs[:,0]],axis=1).max()/DX)
   assert min(dist/DX)>=.1 and relative_contact_step<=g['max_relative_endpoint_step_dx'],'CONTACT_GEOMETRIC_SAFETY'
  Uc,fc,nc=self.contact.evaluate(s['x'],s['d']);self.error('contact_force_N',abs(fc-(s['force']-s['material_force'])).max());self.error('contact_U_J',abs(Uc-s['Uc']))
  self.W+=float(-.5*np.sum((p['force'][bc,2]+s['force'][bc,2])*dx[bc,2])+.5*RHO*V*np.sum(s['velocity'][bc,2]**2-p['velocity'][bc,2]**2))
  dJ=0.
  if len(newfail):
   dist=np.linalg.norm(s['x'][b[newfail,1]]-s['x'][b[newfail,0]],axis=1);dJ=float(.5*self.contact.c*np.sum(np.minimum(dist-DX,0)**2));self.update_components(s)
  self.Jc+=dJ;self.cwork+=float(.5*np.sum((self.prevfc+fc)*dx));D=float(sc[1]);U=float(sc[0]);Know=s['K'];Rnow=self.W-(U-self.initial['scalars'][0]+Know-self.initial['K']+D+s['Uc']-self.initial['Uc']-self.Jc);self.scale=max(self.scale,abs(self.W),abs(U-self.initial['scalars'][0])+abs(Know-self.initial['K'])+D+abs(s['Uc']-self.initial['Uc']));relR=abs(Rnow)/self.scale
  Kfree=.5*RHO*V*np.sum(s['velocity'][mask]**2);dKfree=.5*RHO*V*np.sum(s['velocity'][mask]**2-p['velocity'][mask]**2)
  defect=np.array([dKfree-.5*np.sum((p['force'][mask]+s['force'][mask])*dx[mask]),old['U']-self.last['U']+.5*np.sum((p['material_force']+old['force'])*dx),.5*np.sum((s['material_force']-old['force'])*dx),s['Uc']-p['Uc']-dJ+.5*np.sum((self.prevfc+fc)*dx)]);self.defects+=defect;rem=Rnow+self.defects.sum();self.error('decomposition_J',abs(rem))
  impulse=np.zeros(3);impulse[2]=RHO*V*np.sum(dv[bc,2])-h/2*np.sum(p['force'][bc,2]+s['force'][bc,2]);momentum=float(np.linalg.norm(RHO*V*np.sum(dv,axis=0)-impulse));xeq=float(abs(dx[mask]-h*p['velocity'][mask]-h*h/(2*RHO*V)*p['force'][mask]).max());veq=float(abs(dv[mask]-h/(2*RHO*V)*(p['force'][mask]+s['force'][mask])).max());self.error('momentum_Ns',momentum);self.error('x_VV_m',xeq);self.error('v_VV_m_s',veq)
  u,v=boundary(s['time']);xerr=verr=0.
  for key,sign in [('bottom_ids',-1),('top_ids',1)]:
   ii=self.ids[key]-1;xerr=max(xerr,float(abs(s['x'][ii,2]-self.ref['X'][ii,2]-sign*u).max()));verr=max(verr,float(abs(s['velocity'][ii,2]-sign*v).max()))
  self.error('prescribed_x_m',xerr);self.error('prescribed_v_m_s',verr)
  speed=float(np.linalg.norm(s['velocity'],axis=1).max());disp=float(np.linalg.norm(dx,axis=1).max()/DX);rdisp=float(np.linalg.norm(dx[b[:,1]]-dx[b[:,0]],axis=1).max()/DX);maxdd=float(dd.max());stretch=np.linalg.norm(s['x'][b[:,1]]-s['x'][b[:,0]],axis=1)/self.ref['r0']-1
  assert np.array_equal(s['atom_ids'],np.arange(1,self.N+1)),'ORIGINAL_ATOM_IDENTITY'
  assert np.array_equal(self.ns.field(9),np.arange(1,self.NB+1)),'ORIGINAL_BOND_IDENTITY'
  assert s['native_bonds']==int(np.sum(s['d']<1)) and sc[3]==natural.sum(),'TOPOLOGY_OR_NATURAL_EVENT_COUNT'
  assert sc[4]==s['step'] and sc[5]==s['step'] and abs(s['time']-p['time']-h)<1e-15,'COMMIT_TIME_IDENTITY'
  assert not len(newfail) or np.all(stretch[newfail]>=SC-g['history_abs']),'FAILURE_BELOW_SC'
  assert max(self.maxerr[k] for k in ['kappa','d','history_d_rollback','history_kappa_rollback'])<=g['history_abs'],'IRREVERSIBLE_HISTORY'
  assert sc[2]>=g['no_negative_D_below_J'],'NEGATIVE_D'
  assert max(self.maxerr[k] for k in ['material_force_N','contact_force_N','failed_tension_N'])<=g['force_abs_N'],'FORCE_REFERENCE'
  assert max(self.maxerr[k] for k in ['U_J','D_increment_J','KE_J','contact_U_J','decomposition_J','star_Dsum_J','star_U_sum_J'])<=g['energy_abs_J'],'ENERGY_REFERENCE'
  assert xeq<=g['integrator_equation_abs_m'] and veq<=g['velocity_abs_m_s'] and momentum<=g['momentum_absolute_Ns'],'VERLET_MOMENTUM'
  assert xerr<=g['integrator_equation_abs_m'] and verr<=g['velocity_abs_m_s'],'PRESCRIBED_NORMAL_MOTION'
  assert maxdd<=g['max_damage_increment'],'MAX_DD';assert relR<=g['energy_relative_stop'],'CUMULATIVE_ENERGY_BALANCE'
  assert speed<=g['max_speed_m_s'] and disp<=g['max_step_displacement_dx'] and rdisp<=g['max_relative_endpoint_step_dx'],'MOTION_SAFETY'
  assert np.all(abs(s['x'])<BOX),'DOMAIN_EXIT'
  residual=s['x']-self.Xaug@(self.fit@s['x']);nonaff=float(np.sqrt(np.mean(np.sum(residual**2,axis=1))));Dregion=np.bincount(self.ref['node_regions'],weights=self.Dnode,minlength=len(REGIONS));failed_region=np.bincount(self.ref['bond_regions'][natural],minlength=len(REGIONS));P=RHO*V*s['velocity'].sum(0)
  reaction_bottom=-float(np.sum(s['force'][self.ids['bottom_ids']-1,2]));reaction_top=-float(np.sum(s['force'][self.ids['top_ids']-1,2]));acc=.5*AMPLITUDE*(np.pi/RAMP)**2*np.cos(np.pi*s['time']/RAMP) if s['time']<RAMP else 0.
  row=dict(step=s['step'],time_s=s['time'],dt_s=h,opening_m=2*u,reaction_top_N=reaction_top,reaction_bottom_N=reaction_bottom,generalized_material_reaction_N=(reaction_top-reaction_bottom)/2,actuator_top_N=reaction_top+len(self.ids['top_ids'])*RHO*V*acc,W_J=self.W,U_J=U,K_J=Know,Kfree_J=Kfree,D_J=D,Uc_J=s['Uc'],Jc_J=self.Jc,R_J=Rnow,relative_R=relR,energy_scale_J=self.scale,max_speed_m_s=speed,max_displacement_m=float(np.linalg.norm(s['x']-self.ref['X'],axis=1).max()),nonaffine_RMS_m=nonaff,step_disp_dx=disp,relative_endpoint_step_dx=rdisp,max_dd=maxdd,max_kappa=float(s['kappa'].max()),new_soft=len(newsoft),new_fail=len(newfail),total_natural_failed=int(natural.sum()),initial_cut_count=int(self.ref['initial_cut'].sum()),live_bonds=s['native_bonds'],components=self.components,contact_pairs=nc,max_contact_penetration_dx=penetration,relative_contact_step_dx=relative_contact_step,atom_count=self.N,Px_Ns=P[0],Py_Ns=P[1],Pz_Ns=P[2],momentum_residual_Ns=momentum,x_equation_m=xeq,v_equation_m_s=veq,contact_work_J=self.cwork,VV_defect_J=self.defects[0],spatial_defect_J=self.defects[1],switch_work_J=self.defects[2],contact_defect_J=self.defects[3],remainder_J=rem,min_box_clearance_m=float(np.min(BOX-abs(s['x']))),commits=sc[4],token=sc[5])
  for j,name in enumerate(REGIONS):row['D_'+name+'_J']=Dregion[j];row['failed_'+name]=int(failed_region[j])
  for j,axis in enumerate('xyz'):row[axis+'min_m']=float(s['x'][:,j].min());row[axis+'max_m']=float(s['x'][:,j].max())
  assert np.all(np.isfinite(np.r_[list(row.values()),s['x'].ravel(),s['velocity'].ravel(),s['force'].ravel(),s['kappa'],s['d']])),'NONFINITE'
  self.last=new;self.prevfc=fc;return row,pairs,newsoft,newfail
