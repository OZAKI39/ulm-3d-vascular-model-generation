"""Audit complete original graph cuts of real target/largest free components."""
from common import *
from multifragment_shear_objective import shape

def audit_separation(ref,results):
    graph=module('separation_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
    b=ref['bonds'];cut=ref['initial_cut'];field=load_npz(O/'BIPHASIC_FIELD.npz');rows=[]
    for mode in CASES:
        final_key='final' if results[mode]['native']['status']=='COMPLETED' else 'last_accepted'
        for key in [final_key,'max_multifragment']:
            path=O/'cases'/mode/'keys'/(key+'.npz');s=load_npz(path)
            labels,sizes,attached=graph.topology(ref,s['d']);free_labels=np.flatnonzero(~attached)
            selected=set(np.flatnonzero(~attached&(sizes>=20)&(sizes<=200)).tolist())
            if len(free_labels):selected.add(int(free_labels[np.argmax(sizes[free_labels])]))
            if not selected:
                rows.append(dict(case=mode,key=key,time_s=float(s['time']),status='NO_FREE_COMPONENT_OBSERVED'))
            for label in sorted(selected):
                ids=np.flatnonzero(labels==label);member=np.zeros(len(ref['X']),bool);member[ids]=True
                boundary=member[b[:,0]]^member[b[:,1]];internal=member[b].all(1)
                live_boundary=boundary&(s['d']<1);assert not live_boundary.any()
                natural=boundary&~cut;assert np.all(s['d'][natural]==1)
                internal_live=int((internal&(s['d']<1)).sum());assert internal_live>=len(ids)-1
                rows.append(dict(case=mode,key=key,time_s=float(s['time']),status='PASS_FULL_ORIGINAL_GRAPH_CUT',nodes=len(ids),target_20_200=20<=len(ids)<=200,original_atom_ids=ids+1,
                    all_original_boundary_bonds=int(boundary.sum()),initial_precut_boundary_bonds=int((boundary&cut).sum()),naturally_failed_boundary_bonds=int(natural.sum()),remaining_live_boundary_bonds=int(live_boundary.sum()),
                    naturally_failed_bulk_boundary_bonds=int((natural&(ref['bond_types']==1)).sum()),naturally_failed_interface_boundary_bonds=int((natural&(ref['bond_types']==2)).sum()),naturally_failed_cross_phase_boundary_bonds=int((natural&field['cross_phase']).sum()),
                    internal_live_bonds=internal_live,reference_shape=shape(ref['X'][ids]),current_shape=shape(s['x'][ids]),state_file=str(path.relative_to(O))))
    out=dict(status='PASS',scope='All target-band components and the largest free component at each case last accepted and saved maximum-target state; every original bulk/interface interaction is included.',
        interpretation='Zero live cross-boundary edges certifies discrete material-graph separation. Pre-existing cut edges are counted separately. This does not reconstruct or certify a continuous three-dimensional crack surface.',records=rows)
    dump(O/'GRAPH_SEPARATION_BOUNDARY_AUDIT.json',out)
    return out
