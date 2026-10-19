"""Representative V4 history, contact, mapped-surface and restart checks."""
from common import *
from native import native,predict
from check_elastic import write_patch
from current_surface import CurrentSurface
import time

def material_only(ns):
 ns.l.commands_string(f'group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}')
 ns.dynamic=True;ns.l.command('run 0');return ns.state()

def main():
 start=time.perf_counter();checks={};g=json.loads((O/'ADMISSION_GATES.json').read_text())
 ref=load_npz(O/'tests/elastic/irregular/provenance/reference.npz');case=O/'tests/native_history';write_patch(case,ref)
 ns,C=native(case,install_properties=False,surface=False);model=energy_model(ref)
 # Controlled displacement initializes a small validation coupon only. No V4
 # production topology, thresholds, material identities, or loading are altered.
 ns.static(ref['X']*1.151);p=material_only(ns);s=ns.step(MIN_H)
 expected=p['x']+MIN_H*(p['velocity']+MIN_H*p['force']/(2*ns.mass[:,None]))
 kp,dp=predict(model,ref,expected,p['kappa'],ns.field(13),ns.field(14))
 assert abs(s['x']-expected).max()<2e-18 and abs(s['d']-dp).max()<1e-12
 assert s['d'].max()>0 and (s['d']>=p['d']).all() and (s['kappa']>=p['kappa']).all()
 z=model.evaluate(s['x'],s['d']);before=model.evaluate(s['x'],p['d']);drop=before['U']-z['U']
 assert abs(z['force']-s['material_force']).max()<1e-12 and abs(z['U']-s['native_U'])<1e-17
 assert abs(drop-s['scalars'][1])<1e-17
 direction=z['force']/np.linalg.norm(z['force']);h=1e-9
 fd=(model.evaluate(s['x']+h*direction,s['d'])['U']-model.evaluate(s['x']-h*direction,s['d'])['U'])/(2*h);der=-np.sum(z['force']*direction)
 assert abs(fd-der)/abs(der)<g['energy_gradient_relative_error_max']
 accepted=s;xt=s['x'].copy();xt[:,0]+=1e-7*np.sin(np.arange(ns.N))
 ns.static(xt);restored=ns.static(s['x'])
 for key in ['x','velocity','d','kappa','atom_ids']:assert np.array_equal(restored[key],accepted[key]),key
 assert np.array_equal(restored['scalars'][:6],accepted['scalars'][:6])
 ns.l.command(f'write_restart {case/"accepted.restart"}');np.savez_compressed(case/'accepted.npz',**accepted)
 continuous=ns.step(MIN_H);ns.close()
 rs,_=native(case,install_properties=False,surface=False,restart=case/'accepted.restart');restarted=material_only(rs)
 for key in ['x','velocity','d','kappa','atom_ids']:assert np.array_equal(restarted[key],accepted[key]),key
 r=rs.step(MIN_H);errors={key:float(abs(r[key]-continuous[key]).max()) for key in ['x','velocity','force','d','kappa']}
 errors['physical_scalars']=float(abs(r['scalars'][:6]-continuous['scalars'][:6]).max());rs.close()
 assert errors['x']<1e-18 and errors['velocity']<1e-12 and errors['force']<1e-14 and errors['d']<1e-12 and errors['physical_scalars']<1e-17
 checks['history']=dict(softened_bonds=int((accepted['d']>0).sum()),max_damage=float(accepted['d'].max()),damage_drop_J=drop,softened_force_gradient_relative_error=float(abs(fd-der)/abs(der)),static_trial_rollback_exact=True,native_damage_commits=2,restart_errors=errors)
 # Activate one sub-dx nonuniform-volume contact in the same small coupon.
 cref={k:v.copy() if isinstance(v,np.ndarray) else v for k,v in ref.items()};ib=int(np.argmin(cref['r0']));cref['initial_cut'][ib]=True
 cp=O/'tests/native_contact';write_patch(cp,cref);cn,_=native(cp,install_properties=False,surface=False);cs=cn.static();contact=C(cref);uc,fc,pairs=contact.evaluate(cs['x'],cs['d'])
 assert len(pairs)==1 and uc>0 and abs(uc-cs['Uc'])<1e-17 and abs(cs['force']-cs['material_force']-fc).max()<1e-12
 direction=fc/np.linalg.norm(fc);h=1e-10;fd=(contact.evaluate(cs['x']+h*direction,cs['d'])[0]-contact.evaluate(cs['x']-h*direction,cs['d'])[0])/(2*h);expected=-np.sum(fc*direction)
 assert abs(fd-expected)/abs(expected)<2e-5
 checks['contact']=dict(pairs=pairs+1,nodal_volumes_m3=cref['V'][pairs],native_contact_energy_J=cs['Uc'],independent_energy_J=uc,energy_gradient_relative_error=float(abs(fd-expected)/abs(expected)),net_force_N=fc.sum(0),rule='V3 mean-volume density convention extended symmetrically as ((Vi+Vj)/2)^2; same stiffness and dx onset')
 cn.close()
 # Full-geometry surface audit, affine reproduction and power conservation.
 full=load_npz(O/'provenance/reference.npz');tool=CurrentSurface(full,'uniform');X=full['X'];d=full['initial_cut'].astype(float);v=np.random.default_rng(147).normal(size=X.shape)*.01
 z=tool.evaluate(X,v,d,.001);tool.commit(z);W=z['W'];active=z['active'];force=z['force'];face=z['face_force']
 assert abs(force.sum(0)-face.sum(0)).max()<1e-15
 assert abs(np.cross(X,force).sum(0)-np.cross(z['centroid'],face).sum(0)).max()<1e-18
 assert abs(np.sum(force*v)-np.sum(face*z['face_velocity']))<1e-17 and z['drag_power']<=0
 assert np.max(abs(np.asarray(W.T@X)[active]-z['centroid'][active]))<1e-16
 accepted_top=tool._accepted;ever=tool.ever_loaded.copy();trial=X.copy();trial[:,0]+=.03*X[:,1];tool.geometry(trial,d);tool.rollback_trial()
 assert tool._accepted is accepted_top and np.array_equal(tool.ever_loaded,ever)
 # No physical failure is prescribed: this artificial graph is an isolated
 # mapping unit fixture, used only to exercise the otherwise inactive L2 path.
 mask=(full['regular_X'][:,0]>200e-6)&(full['regular_X'][:,0]<450e-6)&(abs(full['regular_X'][:,1])<125e-6)&(abs(full['regular_X'][:,2])<125e-6)
 b=full['bonds'];fixture=d.copy();fixture[mask[b[:,0]]!=mask[b[:,1]]]=1
 xf=X.copy();xf[mask,0]+=.002;zf=tool.evaluate(xf,v,fixture,.001);l2=zf['active'][tool.L:]
 assert l2.any(),'Mapped L2 fixture must contain supported exposed quadrature'
 face=zf['face_force'];f=zf['force']
 assert abs(f.sum(0)-face.sum(0)).max()<1e-15 and abs(np.cross(xf,f).sum(0)-np.cross(zf['centroid'],face).sum(0)).max()<1e-18
 assert abs(np.sum(f*v)-np.sum(face*zf['face_velocity']))<1e-17
 checks['surface']=dict(reference_active_faces=int(active.sum()),reference_l2_active=int(z['active'][tool.L:].sum()),synthetic_mapping_fixture_l2_active=int(l2.sum()),mapping_fixture_is_not_a_fracture_simulation=True,force_moment_power_conserved=True,accepted_cache_rollback_exact=True,reference_crossing_incidence=tool.bridge.nnz)
 # Real full-geometry native microsteps and one surface-sidecar continuation.
 case=O/'tests/full_surface';ns,_=native(case);p=ns.dynamics();z=ns.tool.evaluate(p['x'],p['velocity'],p['d'],p['time']);ns.tool.commit(z)
 initial=p;Wtot=0.;rows=[]
 for k in range(8):
  s=ns.step(MIN_H);work=ns.tool.work(p,s);Wtot+=work.sum();z=ns.tool.evaluate(s['x'],s['velocity'],s['d'],s['time']);ns.tool.commit(z)
  assert abs(s['external_force']-z['force']).max()<1e-15
  res=Wtot-(s['scalars'][0]+s['K']+s['scalars'][1]+s['Uc']);rows.append(dict(time_s=s['time'],W_J=Wtot,residual_J=res,relative_R=abs(res)/max(abs(Wtot),s['scalars'][0]+s['K']+s['scalars'][1]+s['Uc'],1e-24)))
  if k==6:
   saved=s;ns.l.command(f'write_restart {case/"accepted.restart"}');ns.tool.save_checkpoint(case/'surface.npz',s)
  p=s
 continuous=s;ns.close();rs,_=native(case,restart=case/'accepted.restart');s=rs.dynamics();rs.tool.restore_checkpoint(case/'surface.npz',s)
 for key in ['x','velocity','d','kappa','atom_ids']:assert np.array_equal(s[key],saved[key]),key
 r=rs.step(MIN_H);errors={key:float(abs(r[key]-continuous[key]).max()) for key in ['x','velocity','force','d','kappa']};rs.close()
 assert errors['x']<1e-18 and errors['velocity']<1e-12 and errors['force']<1e-14 and errors['d']<1e-12
 checks['full_surface_native']=dict(microsteps=9,restart_errors=errors,energy_ledger=rows,max_relative_R=max(q['relative_R'] for q in rows))
 checks.update(status='PASS',wall_s=time.perf_counter()-start,library_sha256=sha(LIB));dump(O/'NATIVE_STATE_SURFACE_CHECK.json',checks);print('NATIVE_V4_CHECK_PASS',json.dumps(checks,default=cv),flush=True)
if __name__=='__main__':main()
