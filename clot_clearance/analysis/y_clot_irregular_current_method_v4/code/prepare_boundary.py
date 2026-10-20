from common import *
def main():
    ref=load_npz(O/'provenance/reference.npz');x=ref['X']
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
    print('Y_BOUNDARY_CLASSIFIED',flush=True)
if __name__=='__main__':main()
