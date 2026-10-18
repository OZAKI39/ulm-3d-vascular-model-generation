"""Read-only diagnosis of the seven real V2 end-state member sets."""
from common import *
import csv,time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection,Poly3DCollection

TARGETS={'uniform':[1164,1804,1848,1884],'biphasic':[3769,3809,3846]}

def main():
    started=time.perf_counter();ref=load_npz(V2_O/'provenance/reference.npz')
    old=module('v2_surface_readonly',V2_A/'current_surface.py').CurrentSurface(ref,'uniform')
    b=ref['bonds'];bt=ref['bond_types'];rows=[]
    fig=plt.figure(figsize=(20,16),layout='constrained')
    for mode,cids in TARGETS.items():
        case=V2_O/'cases'/mode;statepath=case/'keys/final.npz'
        s=load_npz(statepath);reg=json.loads((case/'COMPONENT_REGISTRY.json').read_text())
        z=old.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
        assert np.max(abs(z['force']-s['external_force']))<1e-15
        for cid in cids:
            ids=np.array(reg[str(cid)]['atom_ids'],int)-1;member=np.zeros(old.N,bool);member[ids]=True
            internal=member[b].all(1);external=member[b].sum(1)==1;live=s['d']<1
            assert not np.any(external&live)
            face=member[old.tri].all(1);tet=member[old.tet].all(1)
            supported=tet&z['supported_tetra'];candidate=face&z['candidate'];active=face&z['active']
            counts={k:int(np.sum(face&z[k])) for k in ['candidate','original','new_open','accessible','active','unresolved','occluded']}
            counts.update(supported_tetra=int(supported.sum()),all_reference_tetra=int(tet.sum()),degenerate_live_faces=int((face&z['live']&~z['geometric']).sum()),candidate_not_open=int((candidate&~z['new_open']).sum()),opened_not_accessible=int((face&z['new_open']&~z['accessible']).sum()))
            causes=[]
            if not supported.any():causes.append('A_NO_FULLY_LIVE_POSITIVE_VOLUME_AUXILIARY_TETRA')
            if counts['candidate_not_open']:causes.append('B_SUPPORTED_CANDIDATES_FAIL_OPEN_GAP_OR_APERTURE')
            if counts['occluded']:causes.append('C_CURRENT_SOLID_OCCLUSION')
            if counts['opened_not_accessible']:causes.append('D_OPEN_NOT_EXTERIOR_CONNECTED')
            if counts['degenerate_live_faces']:causes.append('E_DEGENERATE_LIVE_TRIANGLES_PRESENT')
            coverage=np.unique(old.tri[active]).size/len(ids)
            rr=dict(case=mode,component_id=cid,nodes=len(ids),original_atom_ids=ids+1,counts=counts,causes=causes,
                active_area_m2=float(z['area'][active].sum()),effective_loaded_area_m2=float((z['area']*z['weight'])[active].sum()),loaded_node_fraction=coverage,
                resultant_force_N=z['force'][ids].sum(0),moment_about_origin_Nm=np.cross(s['x'][ids],z['force'][ids]).sum(0),
                moment_about_COM_Nm=np.cross(s['x'][ids]-s['x'][ids].mean(0),z['force'][ids]).sum(0),nodal_power_W=float(np.sum(z['force'][ids]*s['velocity'][ids])),
                live_bulk_bonds=int((internal&live&(bt==1)).sum()),live_interface_bonds=int((internal&live&(bt==2)).sum()),failed_external_material_bonds=int((external&~live&(bt!=3)).sum()),
                ownership_and_force_mapping='PASS_ORIGINAL_ID_REPLAY; no missing nodal force relative to existing facet formula',
                source_sha256=sha(statepath),source=str(statepath.relative_to(R)),
                candidate_means='Faces adjacent to one surviving positive-volume all-live auxiliary tetra; absent candidates are not evidence of absent material boundary',
                continuum_unloaded_area='UNDETERMINED; node coverage is not continuum wetted-area coverage')
            rows.append(rr)
            selected=face&(z['candidate']|z['original']|z['new_open']|z['active'])
            out=O/'p0_targets'/f'{mode}_{cid}.npz';out.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(out,original_atom_ids=ids+1,X=ref['X'][ids],x=s['x'][ids],velocity=s['velocity'][ids],
                original_live_bulk_bond_ids=np.flatnonzero(internal&live&(bt==1))+1,original_live_interface_bond_ids=np.flatnonzero(internal&live&(bt==2))+1,
                live_bond_atom_ids=b[internal&live]+1,failed_external_bond_ids=np.flatnonzero(external&~live&(bt!=3))+1,
                failed_external_bond_atom_ids=b[external&~live&(bt!=3)]+1,supported_tetra_atom_ids=old.tet[supported]+1,
                face_ids=np.flatnonzero(selected),face_atom_ids=old.tri[selected]+1,**{k:z[k][selected] for k in ['centroid','normal','area','weight','active','original','new_open','accessible','unresolved','aperture_m','occluded','face_velocity','face_force','traction','fluid_region']})
            ax=fig.add_subplot(2,4,len(rows),projection='3d');xyz=s['x']*1e6
            ax.add_collection3d(Line3DCollection(xyz[b[internal&live]],colors='#476d81',linewidths=.5,alpha=.55))
            ax.scatter(*xyz[ids].T,s=11,c='#24576c',depthshade=False)
            for mask,color,alpha in [(candidate&~z['active'],'#e9a23b',.3),(active,'#36ad86',.85)]:
                if mask.any():ax.add_collection3d(Poly3DCollection(xyz[old.tri[mask]],facecolors=color,edgecolors=color,alpha=alpha,linewidths=.8))
            mid=xyz[ids].mean(0);span=np.ptp(xyz[ids],axis=0).max()*1.1
            for j,setlim in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):setlim(mid[j]-span/2,mid[j]+span/2)
            ax.set_box_aspect((1,1,1));ax.view_init(16,-64)
            ax.set_title(f'{mode.title()} ID {cid}: {len(ids)} nodes\n{counts["supported_tetra"]} supported tets; {counts["active"]} active facets\nloaded nodes {100*coverage:.2f}%',fontsize=12)
            ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)');ax.set_zlabel('Z (µm)');ax.tick_params(labelsize=8)
    ax=fig.add_subplot(2,4,8);ax.axis('off');ax.text(0,.95,'Read-only replay of V2 FP64 end states\n\nBlue: real nodes and surviving material bonds\nOrange: supported but unloaded candidates\nGreen: actually loaded V2 facets\n\nOriginal world positions; no inflated shell.\nNo facet does not mean no physical drag.\n\nA/B/C/D/E diagnoses are reported separately.\nForce/ownership replay is audited by original ID.',va='top',fontsize=13,linespacing=1.7)
    fig.suptitle('V2 fragment surface representation: material connectivity versus auxiliary support',fontsize=22)
    (F/'figures').mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(F/'figures'/f'FIG_01_SURFACE_REPRESENTATION_FAILURE.{ext}',dpi=180)
    plt.close(fig)
    result=dict(status='COMPLETED_READ_ONLY_V2_REPLAY',targets=rows,wall_s=time.perf_counter()-started,source_surface_sha256=sha(V2_A/'current_surface.py'))
    dump(O/'P0_TARGET_SURFACE_DIAGNOSIS.json',result);dump(O/'TARGET_COMPONENT_SURFACE_DIAGNOSIS.json',result)
    print(json.dumps([{k:r[k] for k in ['case','component_id','nodes','counts','causes']} for r in rows]),flush=True)

if __name__=='__main__':main()
