"""Distinguish surface vertices from whole fragments lacking represented facets."""
from analyze import *
out={}
ref=load_npz(O/'provenance/reference.npz')
for name in CASES:
    tool=CurrentSurface(ref,name);records=[]
    for path in sorted((O/'cases'/name/'keys').glob('*.npz')):
        if 'failed_unaccepted' in path.name:continue
        s=load_npz(path);z=tool.evaluate(s['x'],s['velocity'],s['d'],float(s['time']));lab,sizes,att=graph.topology(ref,s['d']);tri=tool.tri[z['active']]
        assert np.all(lab[tri]==lab[tri[:,0]][:,None])
        faces=np.bincount(lab[tri[:,0]],minlength=len(sizes));items=[]
        for c in np.flatnonzero(~att):
            ids=np.flatnonzero(lab==c);f=z['force'][ids]
            items.append(dict(original_atom_ids=ids+1,nodes=int(sizes[c]),retained_faces=int(faces[c]),sum_force_norm_N=float(np.linalg.norm(f,axis=1).sum()),sum_drive_force_norm_N=float(np.linalg.norm(z['drive_force'][ids],axis=1).sum()),sum_drag_force_norm_N=float(np.linalg.norm(z['drag_force'][ids],axis=1).sum())))
        records.append(dict(key=path.stem,time_s=float(s['time']),components=items,free_components_with_retained_faces=sum(q['retained_faces']>0 for q in items),free_components_without_retained_faces=sum(q['retained_faces']==0 for q in items),free_nodes_in_components_without_retained_faces=sum(q['nodes'] for q in items if not q['retained_faces'])))
    out[name]=records
dump(O/'FRAGMENT_SURFACE_COVERAGE_AUDIT.json',dict(cases=out,scope='Component-level represented surface coverage, not a complete wetting or hydrodynamic reconstruction. Interior nodes need not be surface vertices. The scalar history free_nodes_without_resolved_surface counts non-vertex nodes, including interior nodes of represented large fragments; it must not be interpreted as that many uncoupled fragments or exposed points.'))
windows=read(O/'MODEL_VALIDITY_WINDOWS.json')
for name,records in out.items():
    final=next(q for q in records if q['key']=='final')
    if 'unresolved_free_surface_nodes_final' in windows[name]:windows[name]['free_nodes_not_vertices_of_retained_facets_final']=windows[name].pop('unresolved_free_surface_nodes_final')
    windows[name]['free_components_without_retained_facets_final']=final['free_components_without_retained_faces']
    windows[name]['free_nodes_in_components_without_retained_facets_final']=final['free_nodes_in_components_without_retained_faces']
dump(O/'MODEL_VALIDITY_WINDOWS.json',windows)
print('COMPONENT_SURFACE_AUDIT', {k:next({a:q[a] for a in ['free_components_with_retained_faces','free_components_without_retained_faces','free_nodes_in_components_without_retained_faces']} for q in v if q['key']=='final') for k,v in out.items()})
