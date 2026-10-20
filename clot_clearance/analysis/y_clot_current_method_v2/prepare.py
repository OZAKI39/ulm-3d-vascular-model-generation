"""Freeze the historical irregular Y geometry and one unoptimized field transfer."""
from common import *
import shutil,datetime
from scipy.ndimage import label,generate_binary_structure

def main():
    assert not (O/'GEOMETRY_FREEZE.json').exists()
    for p in [O/'provenance',O/'build',O/'checks',F]:p.mkdir(parents=True,exist_ok=True)
    for name in ['reference.npz','reference.native','body.data']:
        shutil.copy2(Y_SOURCE/'provenance'/name,O/'provenance'/name)
    ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);x=ref['X'];b=ref['bonds'];bt=ref['bond_types']
    assert n==1927 and len(x)==3488 and len(b)==157637 and not ref['initial_cut'].any()
    for name in ['X','bonds','bond_types','original_bond_ids','r0','taper','m','V','initial_cut']:
        assert np.array_equal(ref[name],load_npz(Y_SOURCE/'provenance/reference.npz')[name]),name
    # One frozen translation, no scaling, phase optimization, random seed search,
    # rebuilding of the Y material graph, or imported rectangular notch.
    translation=np.array([-6.,.05,1.])*1e-3
    old=load_npz(V2_O/'BIPHASIC_FIELD.npz');ix_float=(x[:n]-translation)/DX+np.array([12,12,6]);ix=np.rint(ix_float).astype(int)
    assert np.max(abs(ix_float-ix))<1e-10 and (ix>=0).all() and (ix<np.array(COUNTS)).all()
    domain=old['grid_domain_id'][tuple(ix.T)];strong=domain>0;phase=np.r_[strong.astype(np.int8),np.full(len(x)-n,-1,np.int8)]
    on=S0[bt].copy();off=SC[bt].copy();bulk=bt==1
    node_on=np.where(phase==1,.25,.10);node_off=np.where(phase==1,.30,.12)
    on[bulk]=node_on[b[bulk]].mean(1);off[bulk]=node_off[b[bulk]].mean(1)
    np.savez_compressed(O/'BIPHASIC_FIELD.npz',atom_ids=np.arange(1,len(x)+1),material_phase=phase,
                        strong_domain_id=np.r_[domain,np.zeros(len(x)-n,int)],bond_ids=ref['original_bond_ids'],bond_s0=on,bond_sc=off,
                        translation_m=translation,source_lattice_indices=ix)
    domains=[]
    for k in np.unique(domain[strong]):
        mask=np.zeros(COUNTS,bool);mask[tuple(ix[domain==k].T)]=True
        domains.append(dict(source_domain_id=int(k),nodes=int((domain==k).sum()),face_connected_parts=int(label(mask,generate_binary_structure(3,1))[1])))
    dump(O/'GEOMETRY_FREEZE.json',dict(source=str(Y_SOURCE),material_nodes=n,anchor_nodes=len(x)-n,original_bonds=len(b),
        bulk_bonds=int((bt==1).sum()),interface_bonds=int((bt==2).sum()),scaffold_bonds=int((bt==3).sum()),initial_cut=0,
        geometry_and_reference_graph_exact=True,dx_m=DX,horizon_m=HORIZON,material_bounds_m=[x[:n].min(0),x[:n].max(0)],
        no_rectangular_notch_imported=True,source_sha256={s:sha(Y_SOURCE/'provenance'/s) for s in ['reference.npz','reference.native','body.data']}))
    dump(O/'BIPHASIC_MATERIAL_SPEC.json',dict(label=HET_LABEL,transfer='One rigid translation and geometric clipping of frozen V2 material field',
        translation_m=translation,source_sha256=sha(V2_O/'BIPHASIC_FIELD.npz'),strong_nodes=int(strong.sum()),strong_fraction=float(strong.mean()),
        retained_source_domains=domains,original_rectangular_strong_fraction=2708/8125,scale_changed=False,parameters_searched=False,
        G=G,K=K,rho=RHO,weak=[.10,.12],mixed=[.175,.210],strong=[.25,.30],interface=[S0[2],SC[2]],
        interpretation='Clipping changes phase fraction; this is a transferred application, not an isolated geometry causal comparison.'))
    dump(O/'FROZEN_SHEAR_DRIVE.json',dict(direction=[1,0,0],nominal_Pa=1000,beta_Pa_s_m=BETA,h_eff_m=H_EFF,
        rise_end_s=RISE,hold_end_s=HOLD,drive_off_s=OFF,end_s=END,CFD_forces_applied=False,old_Gaussian_applied=False,
        old_COM_drag_applied=False,numerical_damping_applied=False,acoustic_realizability='UNESTABLISHED'))
    dump(O/'CURRENT_EXPOSED_SURFACE_SPEC.json',dict(L1='Original-grid clipped Freudenthal tetra support, no missing vertices fabricated',
        L2='Frozen current V3 free-component dual faces',lumen='Current centroid and 3um outward probe in original Y tetra lumen',
        wall='Original 20um/50N-per-m/.5 damping ratio; distance-gradient conservative elastic normal; work separately recorded',
        label='PARTIAL_CURRENT_SURFACE_APPROXIMATE_COUPLING',CFD_updated=False))
    # Identify the two genuinely different historical boundaries. A material
    # surface below the vessel wall is not necessarily blood-exposed: the
    # point discretization leaves a geometric offset at the bonded interface.
    from supported_surface import CurrentSurface as SupportedSurface
    from pd_clot.lammps_fragmentation.native_y_wall_query import NativeYWallQuery
    from pd_clot.lammps_fragmentation.y_mural_clot_geometry import locator,hit,make_surface
    from pd_clot.lammps_fragmentation.y_clot_morphology_v3 import profile
    support=SupportedSurface(ref,'uniform');wall=NativeYWallQuery()
    src=R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/V3_GEOMETRY_CONFIG.json'
    config=json.loads(src.read_text());solid,top,attachment=make_surface(wall.sources,config,profile_fn=profile)
    lt=locator(top);la=locator(attachment);kind=np.zeros(len(support.tri),np.int8)
    centers=x[support.tri].mean(1);original=support.original&support.material
    for i in np.flatnonzero(original):
        c=centers[i];end=c+2e-3*support.reference_normal[i];ht=hit(lt,c,end);ha=hit(la,c,end)
        if ht[0] and ha[0]:kind[i]=1 if np.linalg.norm(ht[1]-c)<np.linalg.norm(ha[1]-c) else 2
        elif ht[0]:kind[i]=1
        elif ha[0]:kind[i]=2
    clearance=np.zeros(len(kind));ids=np.flatnonzero(kind==2)
    clearance[ids]=wall.query(centers[ids])['clearance_m']
    np.savez_compressed(O/'BOUNDARY_CLASSIFICATION.npz',tri=support.tri,reference_class=kind,reference_wall_clearance_m=clearance)
    dump(O/'BOUNDARY_CLASSIFICATION.json',dict(exposed_faces=int((kind==1).sum()),attachment_faces=int((kind==2).sum()),
        unresolved_original_faces=int((original&(kind==0)).sum()),source_config_sha256=sha(src),
        reconstructed_exposed_area_m2=float(top.area),reconstructed_attachment_area_m2=float(attachment.area),
        policy='Original attachment facets never load while anchor-connected. A free component requires an additional current wall gap >1um, smooth to5um, and ordinary solid/lumen checks.'))
    print('Y_GEOMETRY_FROZEN',n,len(x)-n,int(strong.sum()),float(strong.mean()),domains,flush=True)

if __name__=='__main__':main()

