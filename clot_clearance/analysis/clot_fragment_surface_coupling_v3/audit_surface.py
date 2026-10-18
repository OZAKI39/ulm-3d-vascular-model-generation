"""Small read-only audits of the new surface quadrature on real V2 states."""
from common import *
from current_surface import CurrentSurface,CLASS_NAMES
from diagnose_v2 import TARGETS
import time,csv
from scipy.spatial.transform import Rotation

def main():
    start=time.perf_counter();ref=load_npz(V2_O/'provenance/reference.npz');surf=CurrentSurface(ref,'uniform');rows=[];checks=[]
    for mode,cids in TARGETS.items():
        case=V2_O/'cases'/mode;s=load_npz(case/'keys/final.npz');reg=json.loads((case/'COMPONENT_REGISTRY.json').read_text())
        z=surf.evaluate(s['x'],s['velocity'],s['d'],float(s['time']));old=surf.base.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
        active=z['active'];W=z['W'];mapped=W.T@s['x'];ff=z['face_force']
        err=dict(centroid_m=float(np.max(abs(mapped[active]-z['centroid'][active]),initial=0)),force_N=float(np.max(abs(z['force'].sum(0)-ff.sum(0)))),
                 moment_Nm=float(np.max(abs(np.cross(s['x'],z['force']).sum(0)-np.cross(z['centroid'],ff).sum(0)))),power_W=float(abs(np.sum(s['velocity']*z['force'])-np.sum(z['face_velocity']*ff))))
        assert err['centroid_m']<1e-15 and err['force_N']<1e-15 and err['moment_Nm']<1e-17 and err['power_W']<1e-15,err
        counts={k:int(np.sum((z['face_class'][surf.L:]==i)&z['candidate'][surf.L:])) for i,k in enumerate(CLASS_NAMES)}
        checks.append(dict(case=mode,mapping=err,classes_L2=counts,topology_builds=surf.topology_builds))
        w2=W[:,surf.L:][:,active[surf.L:]]
        checks[-1]['affine_weights']=dict(minimum=float(w2.data.min(initial=0)),maximum_absolute_column_sum=float(np.asarray(abs(w2).sum(0)).max(initial=0)),negative_entries=int((w2.data<0).sum()),total_nonzero_entries=int(w2.nnz))
        for cid in cids:
            ids=np.array(reg[str(cid)]['atom_ids'])-1;member=np.zeros(surf.N,bool);member[ids]=True
            oldmask=member[surf.base.tri].all(1);mask=np.r_[oldmask,member[surf.owner]]
            loaded=mask&active;cols=np.flatnonzero(loaded);covered=np.unique(W[:,cols].nonzero()[0]);covered=covered[member[covered]]
            candidates=mask&z['candidate'];unknown=mask&z['unresolved'];l2=np.arange(len(mask))>=surf.L
            rr=dict(case=mode,component_id=cid,nodes=len(ids),old_area_m2=float(np.sum(old['area'][oldmask]*old['weight'][oldmask])),new_area_m2=float(np.sum(z['area'][mask]*z['weight'][mask])),
                original_exposed_area_m2=float(np.sum(z['area'][mask&z['original']]*z['weight'][mask&z['original']])),new_crack_area_m2=float(np.sum(z['area'][mask&z['accessible']]*z['weight'][mask&z['accessible']])),
                L1_active_facets=int((loaded&~l2).sum()),L2_active_patches=int((loaded&l2).sum()),candidate_patches=int(candidates.sum()),unresolved_patches=int(unknown.sum()),
                unresolved_candidate_fraction=float(unknown.sum()/max(1,candidates.sum())),old_node_coverage=np.unique(surf.base.tri[oldmask&old['active']]).size/len(ids),new_node_coverage=len(covered)/len(ids),
                L2_reference_bridge_exclusions=int((member[surf.owner]&z['topology']['boundary']&(z['l2_bridge_count']>0)).sum()),
                L2_existing_facet_overlap_exclusions=int((member[surf.owner]&z['topology']['boundary']&z['l2_overlap_guard']).sum()),
                force_norm_N=float(np.linalg.norm(z['force'][ids].sum(0))),moment_COM_Nm=float(np.linalg.norm(np.cross(s['x'][ids]-s['x'][ids].mean(0),z['force'][ids]).sum(0))),power_W=float(np.sum(s['velocity'][ids]*z['force'][ids])))
            rows.append(rr)
            np.savez_compressed(O/'p0_targets'/f'{mode}_{cid}_v3_quadrature.npz',original_atom_ids=ids+1,face_ids=np.flatnonzero(mask),
                **{k:z[k][mask] for k in ['centroid','normal','area','active','weight','face_velocity','outer_velocity','face_force','traction','face_class','candidate','unresolved']},method=surf.method[mask])
        # Exact rigid translation and rigid spin compatibility, on fixed real geometry.
        u=np.array([.3,-.2,.1]);v=np.broadcast_to(u,s['velocity'].shape).copy()
        t=END;zz=surf.evaluate(s['x'],v,s['d'],t)
        projected=u-z['normal']*np.sum(u*z['normal'],axis=1)[:,None]
        expected=-BETA*(z['area']*z['weight'])[:,None]*projected
        translation=float(np.max(abs(zz['face_force']-expected)));assert translation<1e-15
        uniform=surf.evaluate(s['x'],np.zeros_like(v),s['d'],HOLD)
        ex=np.array([1.,0.,0.]);tangent=ex-z['normal']*z['normal'][:,0,None]
        expected_shear=1000*(z['area']*z['weight'])[:,None]*tangent
        constant_shear=float(np.max(abs(uniform['force'].sum(0)-expected_shear.sum(0))));assert constant_shear<1e-14
        omega=np.array([23.,-17.,31.]);vr=np.cross(omega,s['x']);zz=surf.evaluate(s['x'],vr,s['d'],t)
        rotation_velocity=float(np.max(abs(zz['face_velocity'][active]-np.cross(omega,z['centroid'][active]))));assert rotation_velocity<1e-12
        Q=Rotation.from_rotvec([.21,-.34,.45]).as_matrix();shift=np.array([.001,-.002,.003]);xx=s['x']@Q.T+shift
        rotated=surf.geometry(xx,s['d']);common=active&rotated['active'];changed=int(np.sum(active!=rotated['active']))
        area_error=float(np.max(abs(rotated['area'][common]-z['area'][common])));normal_error=float(np.max(abs(rotated['normal'][common]-z['normal'][common]@Q.T)))
        assert area_error<1e-17 and normal_error<1e-10
        # L1 triangle ray/probe proxy has its own pre-existing roundoff edge cases.
        checks[-1].update(constant_1000Pa_shear_resultant_error_N=constant_shear,rigid_translation_force_error_N=translation,rigid_spin_velocity_error_m_s=rotation_velocity,
                         rotated_area_error_m2=area_error,rotated_normal_error=normal_error,rotated_active_face_disagreement_count=changed)
        surf.commit(z);accepted=surf._accepted;ever=surf.ever_loaded.copy();total=surf.cumulative_new_area
        trial=s['d'].copy();trial[np.flatnonzero(trial<1)[:2]]=1;surf.geometry(s['x'],trial);surf.rollback_trial()
        assert surf._accepted is accepted and np.array_equal(ever,surf.ever_loaded) and total==surf.cumulative_new_area
        checks[-1]['rejected_surface_trial_did_not_commit']=True
    with (O/'SURFACE_QUADRATURE_AUDIT.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    dump(O/'SURFACE_GEOMETRY_AND_MAPPING_AUDIT.json',dict(status='PASS',checks=checks,targets=rows,wall_s=time.perf_counter()-start,
        limitations=['Coverage uses actually nonzero quadrature weights, not true continuum area','L2 ray-accessible affine dual cells are an approximate partial boundary','Bounded signed affine weights preserve moment/power; nonnegative extrapolation is not forced','Rank-deficient/bridged/overlapping/occluded cells stay unresolved']))
    print(json.dumps(dict(checks=checks,targets=rows,wall_s=time.perf_counter()-start)),flush=True)

if __name__=='__main__':main()
