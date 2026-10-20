"""Read an immutable accepted stage from both cases without advancing dynamics."""
from common import *
from current_surface import CurrentSurface
from multifragment_shear_objective import metrics, shape
import argparse

def main(argv=None):
    p=argparse.ArgumentParser();selection=p.add_mutually_exclusive_group(required=True);selection.add_argument('--key',choices=['peak_load','upper_drive_off','end_state']);selection.add_argument('--time-ms',type=float);p.add_argument('--nearest',action='store_true',help='Use the nearest real frame within half of the 25 us saving interval; record its actual time.');args=p.parse_args(argv)
    assert not args.nearest or args.time_ms is not None
    tag=args.key if args.key else (f'time_{args.time_ms:.3f}_ms'+('_nearest' if args.nearest else '')).replace('.','p')
    ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz');b=ref['bonds'];n=int(ref['clot_count'])
    graph=module('stage_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py');out={}
    for mode in CASES:
        keys=O/'cases'/mode/'keys'
        summary_path=keys.parent/'SUMMARY.json';summary=json.loads(summary_path.read_text()) if summary_path.exists() else {}
        registry_path=keys.parent/'COMPONENT_REGISTRY.json'
        registry_lookup={tuple(r['atom_ids']):int(cid) for cid,r in json.loads(registry_path.read_text()).items()} if registry_path.exists() else {}
        if args.key:
            actual_key=args.key
            if args.key=='end_state':
                assert summary,'End-state replay requires the finished native summary'
                actual_key='final' if summary['status']=='COMPLETED' else 'last_accepted'
            path=keys/(actual_key+'.npz')
        else:
            matches=[];available=[]
            for candidate in keys.glob('frame_*.npz'):
                if candidate.name.endswith('.tmp.npz'):continue
                with np.load(candidate,allow_pickle=False) as data:
                    gap=abs(float(data['time'])*1e3-args.time_ms);available.append((gap,candidate))
                    if gap<1e-8:matches.append(candidate)
            if args.nearest:
                gap,candidate=min(available)
                assert gap<=.0125+1e-8,('NO_NEARBY_ACCEPTED_FRAME',mode,gap)
                matches=[candidate]
            assert len(matches)==1,('IMMUTABLE_ACCEPTED_FRAME_NOT_OBSERVED',mode,args.time_ms,matches)
            path=matches[0]
        s=load_npz(path);lab,sizes,attached=graph.topology(ref,s['d']);free=np.flatnonzero(~attached)
        surf=CurrentSurface(ref,'uniform');z=surf.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
        assert np.max(np.abs(z['force']-s['external_force']))<1e-15
        records=[]
        selected=sorted(free,key=lambda i:-sizes[i])[:5]
        selected=sorted(set(selected)|set(np.flatnonzero(~attached&(sizes>=20))),key=lambda i:-sizes[i])
        for label in selected:
            ids=np.flatnonzero(lab==label);member=np.zeros(len(ref['X']),bool);member[ids]=True
            internal=member[b].all(1)&(s['d']<1);boundary=member[b[:,0]]^member[b[:,1]]
            assert not (boundary&(s['d']<1)).any()
            faces=z['active']&member[surf.tri].all(1);domains=field['strong_domain_id'][ids]
            records.append(dict(component_id=registry_lookup.get(tuple((ids+1).tolist())),nodes=len(ids),original_atom_ids=ids+1,internal_alive_bonds=int(internal.sum()),live_outgoing_bonds=0,strong_nodes=int((field['material_phase'][ids]==1).sum()),strong_domain_node_counts={int(c):int((domains==c).sum()) for c in np.unique(domains) if c>0},live_internal_cross_phase_bonds=int((internal&field['cross_phase']).sum()),live_internal_interdomain_strong_bonds=int((internal&(field['material_phase'][b]==1).all(1)&(field['strong_domain_id'][b[:,0]]!=field['strong_domain_id'][b[:,1]])).sum()),reference_shape=shape(ref['X'][ids]),actual_shape=shape(s['x'][ids]),active_facets=int(faces.sum()),active_area_m2=float(z['area'][faces].sum()),nodes_on_active_faces=len(np.unique(surf.tri[faces])),actual_resultant_force_N=s['external_force'][ids].sum(0)))
        if args.key=='end_state':assert all(r['component_id'] is not None for r in records),'End-state members must resolve to the native original-ID registry'
        out[mode]=dict(key=tag,native_status=summary.get('status','NOT_YET_SUMMARIZED'),time_s=float(s['time']),requested_time_ms=args.time_ms,time_alignment='NEAREST_REAL_ACCEPTED_FRAME_NO_INTERPOLATION' if args.nearest else 'EXACT_STAGE_OR_FRAME',step=int(s['step']),source_path=str(path.relative_to(O)),source_sha256=sha(path),metrics=metrics(sizes[free],n),largest_free_sizes=sorted(sizes[free].tolist(),reverse=True)[:15],components=records)
    out['scope']='Immutable actual accepted stage comparison, not a final-window conclusion. All original live material interactions determine connectivity; U phase labels are spatial covariates only.'
    dump(O/(tag.upper()+'_CASE_COMPARISON_READBACK.json'),out)
    print(json.dumps({m:dict(t_ms=out[m]['time_s']*1e3,largest=out[m]['largest_free_sizes'],targets=out[m]['metrics']['n_fragments_20_200'],largest_component={k:v for k,v in out[m]['components'][0].items() if k in ['nodes','internal_alive_bonds','strong_nodes','strong_domain_node_counts','live_internal_cross_phase_bonds','live_internal_interdomain_strong_bonds','active_facets']}) for m in CASES},default=lambda x:x.tolist()),flush=True)

if __name__=='__main__':main()
