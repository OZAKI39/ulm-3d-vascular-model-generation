from native import *
from typed_model import TypedLPS
from external import CONTACT
class Mechanics:
 def __init__(self,ns,initial,gates):
  self.ns=ns;self.ref=ns.ref;self.b=self.ref['bonds'];self.bt=self.ref['bond_types'];self.n=ns.n;self.N=ns.N;self.g=gates;self.initial=initial;self.model=TypedLPS(self.ref);self.contact=Contact(self.ref);self.last=self.model.evaluate(initial['x'],initial['d']);self.prevfc=self.contact.evaluate(initial['x'],initial['d'])[1];self.work={k:0. for k in ['load','drag','numeric','wall_elastic','wall_damping']};self.Jc=0.;self.scale=1e-15;self.defects=np.zeros(5);self.maxerr={};self.flow_supply=0.;self.drag_dissipation=0.;self.last_row=None
  self.Uw0=.5*CONTACT['stiffness_N_m']*np.sum(initial['wall_overlap']**2)
 def error(self,k,v):self.maxerr[k]=max(self.maxerr.get(k,0.),float(v))
 def check(self,p,s,h):
  b=self.b;bt=self.bt;n=self.n;g=self.g;dx=s['x']-p['x'];dv=s['velocity']-p['velocity'];dd=s['d']-p['d'];sc=s['scalars'];newfail=np.flatnonzero((s['d']==1)&(p['d']<1));newsoft=np.flatnonzero((s['d']>0)&(p['d']==0))
  new=self.model.evaluate(s['x'],s['d']);old=self.model.evaluate(s['x'],p['d']) if np.any(dd) else new;kk,dr=self.model.update(s['x'],p['kappa'])
  self.error('material_force_N',abs(new['force']-s['material_force']).max());self.error('U_J',abs(new['U']-sc[0]));self.error('D_increment_J',abs(old['U']-new['U']-(sc[1]-p['scalars'][1])));self.error('kappa',abs(kk-s['kappa']).max());self.error('d',abs(dr-s['d']).max());self.error('history_d_rollback',np.max(-dd));self.error('history_kappa_rollback',np.max(p['kappa']-s['kappa']));self.error('failed_tension_N',s['failed_tension_max_N'])
  self.error('KE_J',abs(.5*MASS*np.sum(s['velocity'][:n]**2)-s['K']))
  Uc,fc,pairs=self.contact.evaluate(s['x'],s['d']);external=np.zeros_like(fc);external[:n]=s['external_force'];self.error('contact_force_N',abs(s['force']-s['material_force']-external-fc).max());self.error('contact_U_J',abs(Uc-s['Uc']))
  for key in self.work:self.work[key]+=float(.5*np.sum((p[key]+s[key])*dx[:n]))
  # Relative drag loss and background supply are diagnostic half-velocity time
  # integrals. The primary mechanical ledger uses actual force-displacement work.
  for z in [p,s]:
   self.flow_supply+=h/2*float(np.sum(z['drag']*z['fluid_velocity']))
   self.drag_dissipation-=h/2*float(np.sum(z['drag']*(z['velocity_force_half']-z['fluid_velocity'])))
  dJ=0.
  if len(newfail):
   rr=np.linalg.norm(s['x'][b[newfail,1]]-s['x'][b[newfail,0]],axis=1);dJ=float(.5*self.contact.c*np.sum(np.minimum(rr-DX,0)**2))
  self.Jc+=dJ;D=float(sc[1]);U=float(sc[0]);Know=s['K'];Uw=float(.5*CONTACT['stiffness_N_m']*np.sum(s['wall_overlap']**2));Up=float(.5*CONTACT['stiffness_N_m']*np.sum(p['wall_overlap']**2))
  W=sum(self.work[k] for k in ['load','drag','numeric','wall_damping']);deltaE=U-self.initial['scalars'][0]+Know-self.initial['K']+D+s['Uc']-self.initial['Uc']+Uw-self.Uw0-self.Jc;Rnow=W-deltaE
  self.scale=max(self.scale,abs(W),sum(abs(v) for k,v in self.work.items() if k!='wall_elastic'),abs(U-self.initial['scalars'][0])+abs(Know-self.initial['K'])+D+abs(s['Uc']-self.initial['Uc'])+abs(Uw-self.Uw0));rel=abs(Rnow)/self.scale
  defect=np.array([s['K']-p['K']-.5*np.sum((p['force'][:n]+s['force'][:n])*dx[:n]),old['U']-self.last['U']+.5*np.sum((p['material_force']+old['force'])*dx),.5*np.sum((s['material_force']-old['force'])*dx),s['Uc']-p['Uc']-dJ+.5*np.sum((self.prevfc+fc)*dx),Uw-Up+.5*np.sum((p['wall_elastic']+s['wall_elastic'])*dx[:n])]);self.defects+=defect;self.error('decomposition_J',abs(Rnow+self.defects.sum()))
  xeq=float(abs(dx[:n]-h*p['velocity'][:n]-h*h/(2*MASS)*p['force'][:n]).max());veq=float(abs(dv[:n]-h/(2*MASS)*(p['force'][:n]+s['force'][:n])).max());self.error('x_VV_m',xeq);self.error('v_VV_m_s',veq)
  reactions=-s['force'][n:].sum(0);oldreaction=-p['force'][n:].sum(0);impulse=h/2*(s['external_force'].sum(0)+p['external_force'].sum(0)+reactions+oldreaction);momentum=float(np.linalg.norm(MASS*np.sum(dv[:n],axis=0)-impulse));self.error('momentum_Ns',momentum)
  speed=float(np.linalg.norm(s['velocity'][:n],axis=1).max());stepdisp=float(np.linalg.norm(dx[:n],axis=1).max()/DX);rdisp=float(np.linalg.norm(dx[b[:,1]]-dx[b[:,0]],axis=1).max()/DX);maxdd=float(dd.max());pen=0.;contactstep=0.
  if len(pairs):
   rr=np.linalg.norm(s['x'][pairs[:,1]]-s['x'][pairs[:,0]],axis=1);pen=float(np.max(1-rr/DX));contactstep=float(np.linalg.norm(dx[pairs[:,1]]-dx[pairs[:,0]],axis=1).max()/DX)
   assert min(rr/DX)>=.1 and contactstep<=g['max_relative_endpoint_step_dx'],'CONTACT_GEOMETRIC_SAFETY'
  det=s['detached'];sz=s['component_sizes'];att=s['component_attached'];free=sz[~att];largest=int(free.max(initial=0));tiny=int(free[free<=2].sum())
  assert int(sz.sum())==n and int(free.sum())==int(det.sum()),'MATERIAL_MASS_CONSERVATION'
  assert np.array_equal(s['atom_ids'],np.arange(1,self.N+1)) and np.array_equal(self.ns.field(9),np.arange(1,len(b)+1)) and np.array_equal(self.ns.field(12),bt),'ORIGINAL_IDENTITIES'
  assert np.array_equal(s['x'][n:],self.ref['X'][n:]) and not s['velocity'][n:].any(),'FIXED_ANCHORS'
  assert not s['d'][bt==3].any(),'FIXED_ANCHOR_SCAFFOLD_DAMAGE'
  assert not np.any(s['drag'][~det]) and not np.any(s['numeric'][det]),'ATTACHMENT_FORCE_SEMANTICS'
  assert not np.any(p['detached']&~det),'IRREVERSIBLE_DETACHMENT'
  assert s['native_bonds']==int(np.sum(s['d']<1)) and sc[3]==int(np.sum(s['d']==1)),'TOPOLOGY_COUNT'
  if h:assert s['step']==p['step']+1 and sc[4]==s['step'] and sc[5]==s['step'] and abs(s['time']-p['time']-h)<1e-15,'DYNAMIC_COMMIT_TIME'
  assert max(self.maxerr[k] for k in ['kappa','d','history_d_rollback','history_kappa_rollback'])<=g['history_abs'],'IRREVERSIBLE_HISTORY'
  assert sc[2]>=g['no_negative_D_below_J'],'NEGATIVE_D'
  assert max(self.maxerr[k] for k in ['material_force_N','contact_force_N','failed_tension_N'])<=g['force_abs_N'],'FORCE_REFERENCE'
  assert max(self.maxerr[k] for k in ['U_J','D_increment_J','KE_J','contact_U_J','decomposition_J'])<=g['energy_abs_J'],'ENERGY_REFERENCE'
  assert xeq<=g['integrator_equation_abs_m'] and veq<=g['velocity_abs_m_s'] and momentum<=g['momentum_absolute_Ns'],'VERLET_MOMENTUM'
  assert maxdd<=g['max_damage_increment'] and rel<=g['energy_relative_stop'],('DAMAGE_OR_ENERGY_GATE',maxdd,rel)
  assert speed<=g['max_speed_m_s'] and stepdisp<=g['max_step_displacement_dx'] and rdisp<=g['max_relative_endpoint_step_dx'],'MOTION_SAFETY'
  row=dict(step=s['step'],time_s=s['time'],dt_s=h,load_amp=self.ns.ext.amp(s['time']),load_sum_norm_N=float(np.linalg.norm(s['load'],axis=1).sum()),load_max_N=float(np.linalg.norm(s['load'],axis=1).max()),max_speed_m_s=speed,max_step_displacement_dx=stepdisp,max_relative_step_dx=rdisp,max_dd=maxdd,max_kappa=float(s['kappa'].max()),new_soft=len(newsoft),new_fail=len(newfail),bulk_softened=int(np.sum((bt==1)&(s['d']>0))),interface_softened=int(np.sum((bt==2)&(s['d']>0))),bulk_failed=int(np.sum((bt==1)&(s['d']==1))),interface_failed=int(np.sum((bt==2)&(s['d']==1))),live_bonds=s['native_bonds'],attached_nodes=int(n-det.sum()),detached_nodes=int(det.sum()),attached_mass_kg=float((n-det.sum())*MASS),detached_mass_kg=float(det.sum()*MASS),detached_components=len(free),largest_detached_nodes=largest,largest_detached_diameter_m=float((6*largest*V/np.pi)**(1/3)),tiny_components=int(np.sum(free<=2)),tiny_nodes=tiny,unresolved_debris_mass_kg=tiny*MASS,N_ge200=int(np.sum(free*V>=np.pi*(200e-6)**3/6)),N_ge300=int(np.sum(free*V>=np.pi*(300e-6)**3/6)),valid_drag_nodes=int(s['valid_flow'].sum()),invalid_detached_nodes=int(np.sum(det&~s['valid_flow'])),W_J=W,U_J=U,K_J=Know,D_J=D,Uc_J=s['Uc'],Uw_J=Uw,Jc_J=self.Jc,R_J=Rnow,relative_R=rel,energy_scale_J=self.scale,drag_background_supply_J=self.flow_supply,drag_relative_loss_J=self.drag_dissipation,contact_pairs=len(pairs),max_contact_penetration_dx=pen,max_wall_overlap_m=float(s['wall_overlap'].max()),min_wall_clearance_m=float(s['clearance'].min()),wall_contact_nodes=int(np.sum(s['wall_overlap']>0)),exited_nodes=int(s['exited'].sum()),momentum_residual_Ns=momentum,x_equation_m=xeq,v_equation_m_s=veq,atom_count=self.N,commits=sc[4],token=sc[5])
  for key,value in self.work.items():row['W_'+key+'_J']=value
  for i,key in enumerate(['VV','spatial','switch','particle_contact','wall_contact']):row[key+'_defect_J']=self.defects[i]
  for i,axis in enumerate('xyz'):row['anchor_reaction_'+axis+'_N']=float(reactions[i]);row['load_'+axis+'_N']=float(s['load'][:,i].sum())
  assert np.all(np.isfinite(np.r_[list(row.values()),s['x'].ravel(),s['velocity'].ravel(),s['force'].ravel(),s['kappa'],s['d']])),'NONFINITE'
  self.last=new;self.prevfc=fc;self.last_row=row;return row,pairs,newsoft,newfail
