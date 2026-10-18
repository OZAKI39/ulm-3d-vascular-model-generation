"""Distinguish physical quadrature sites from mapped nodal support coverage."""
from common import *
from current_surface import CurrentSurface
from diagnose_v2 import TARGETS
from analyze import write_csv

def main():
    ref=load_npz(V2_O/'provenance/reference.npz');surf=CurrentSurface(ref,'uniform');rows=[]
    for mode,cids in TARGETS.items():
        case=V2_O/'cases'/mode;s=load_npz(case/'keys/final.npz');reg=json.loads((case/'COMPONENT_REGISTRY.json').read_text())
        z=surf.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
        for cid in cids:
            ids=np.array(reg[str(cid)]['atom_ids'])-1;member=np.zeros(surf.N,bool);member[ids]=True
            mask=member[surf.tri[:,0]]&z['active'];l1=mask[:surf.L];l2=mask[surf.L:]
            sites=np.unique(np.r_[surf.base.tri[l1].ravel(),surf.owner[l2]])
            mapped=np.unique(z['W'][:,mask].nonzero()[0]);assert np.all(member[mapped])
            rows.append(dict(case=mode,component_id=cid,nodes=len(ids),L1_active_facets=int(l1.sum()),L2_active_dual_faces=int(l2.sum()),
                distinct_quadrature_owner_or_vertex_nodes=len(sites),quadrature_owner_or_vertex_fraction=len(sites)/len(ids),
                mapped_force_support_nodes=len(mapped),mapped_force_support_fraction=len(mapped)/len(ids),
                effective_loaded_area_m2=float(np.sum(z['area'][mask]*z['weight'][mask])),true_continuum_surface_coverage='UNKNOWN'))
    write_csv(O/'SURFACE_COVERAGE_DEFINITIONS.csv',rows)
    dump(O/'SURFACE_COVERAGE_DEFINITIONS.json',dict(status='AUDITED',targets=rows,
        quadrature_site_definition='Distinct L1 vertices plus owners of L2 dual faces; still not a continuum surface-area fraction',
        mapping_support_definition='Distinct original nodes with nonzero shape weights on active quadrature elements; may include interior support nodes',
        warning='A broader force distribution is not proof of a more complete physical shell; both metrics and unresolved candidates must be retained'))
    print('COVERAGE_DEFINITIONS_COMPLETE',rows,flush=True)

if __name__=='__main__':main()
