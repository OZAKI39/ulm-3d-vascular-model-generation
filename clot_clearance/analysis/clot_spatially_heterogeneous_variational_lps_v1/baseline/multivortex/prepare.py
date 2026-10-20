from common import *
from current_surface import CurrentSurface,OPEN,FULL_OPEN
from multifragment_shear_objective import freeze
from work_controller import SPEC
import shutil,subprocess,sys,datetime

assert not (O/'BASELINE_MANIFEST.json').exists()
assert sha(LIB)=='465f8144c6da5c9de3aec4137b9c325f53c26737b199c536fb4a20845511ef77'
(O/'provenance').mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
for name in ['reference.npz','reference.native','body.data']:shutil.copy2(OLD/'provenance'/name,O/'provenance'/name)
ref=load_npz(O/'provenance/reference.npz');surf=CurrentSurface(ref,'counter_rotating')
np.savez_compressed(O/'provenance/auxiliary_surface_geometry.npz',tet_atom_ids=surf.tet+1,triangle_atom_ids=surf.tri+1,
    face_owners=surf.owners,face_opposite_atom_ids=surf.opposite+1,tetra_bond_ids=surf.te+1,
    face_edge_bond_ids=surf.edge_ids+1,original=surf.original,material_faces=surf.material,
    reference_area_m2=surf.reference_area,reference_normal=surf.reference_normal)
deps=[LIB,OLD/'DELIVERY_MANIFEST.json',R/'analysis/v3_softening_pilot_v1/native.py',
    R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',PA/'native.py']
dump(O/'BASELINE_MANIFEST.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files={str(p):dict(sha256=sha(p),bytes=p.stat().st_size) for p in deps}))
shutil.copy2(OLD/'MATERIAL_AND_SPECIMEN_FREEZE.json',O/'MATERIAL_AND_SPECIMEN_FREEZE.json')
dump(O/'CURRENT_EXPOSED_SURFACE_SPEC.json',dict(status='IMPLEMENTED_CONSERVATIVE_GEOMETRIC_PROXY',label='QUASISTATIC_FLUID_ACCESSIBILITY_PROXY',
    tetra_count=len(surf.tet),unique_triangle_count=len(surf.tri),method='Freudenthal six-tetra subdivision of the actual reference lattice, including fixed anchors as occlusion geometry',
    support='All six original tetra bonds live, positive signed volume ratio; surface has three actual original material IDs and live edge bonds. Auxiliary tetrahedra carry no mass/stiffness and do not change native mechanics.',
    gap='Six-edge local graph disconnected AND positive volume AND (Vcurrent/Vref-1)*dx > aperture_min. Failed support alone is not fluid.',
    portal='Shared actual triangular face local graph disconnected AND area-dilation proxy (A/Aref-1)*dx > aperture_min; gap cells connect through such portals to an outer portal.',
    new_face='Boundary of a fully supported material tetra against an exterior-connected open gap; opposing node normal-height increase > aperture_min. Each unique triangle mapped once.',
    aperture_min_m=OPEN,traction_full_aperture_m=FULL_OPEN,regularization='C1 smooth aperture gate from 1 to 5 um, same factor on drive and drag',
    categories=['ORIGINAL_EXPOSED_SURFACE','NEW_OPEN_CRACK_SURFACE','FLUID_ACCESSIBLE_CRACK_SURFACE','DETACHED_FRAGMENT_EXPOSED_SURFACE','CLOSED_OR_UNRESOLVED_CRACK_SURFACE'],
    classification='original/new/open/accessibility are geometry flags; detached is orthogonal actual native component flag, not a second added force',
    geometric_limitations=['Node-centered support-complex boundary is a resolution-dependent proxy, not unique continuum crack location.','Damaged connected cells, inverted/closed cells, sheet/chain/fine components without 3D supported tetra remain unresolved.','Local void transmission ignores nonlocal bond lines as impermeable solid tubes; long live bridges still determine material connectivity. No fiber radii or pore flow model supplied.','No guaranteed lower bound on true continuum wetted area; unresolved true area cannot be measured from point data alone. Report candidate-area and component coverage, never claim zero unknown area.','Current gap connectivity uses reference adjacency; remote self-contact can occlude it and is not resolved by this proxy. Must audit closure/overlap before admission.'],
    mass='Native original node IDs, volumes and mass are untouched; auxiliary volume is NOT native material volume.',
    transaction='Geometry is pure current-state recomputation; only accepted commit changes cumulative exposure. No rejected predictor changes damage/surface history.'))
dump(O/'FRAGMENT_FLUID_TRACTION.json',dict(law='tau=(mu/h_eff)*P_t*(u_outer-v_face)',mu=MU,h_eff=H_EFF,mapping='positive 1/3 to original triangle vertices',
    detached='same spatial and relative-velocity law on supported reachable faces; no attachment force switch',
    fines='No artificial closed shell and no point drag added; unrepresented hydrodynamics remains unresolved',
    scope='Hypothetical surface streamfunction field on each oriented facet. Not a single solved global 3D liquid velocity or acoustic streaming. Opposite face normals affect tangent circulation.'))
dump(O/'WORK_MATCH_SPEC.json',SPEC);freeze()
cells=[];lo=ref['X'][:int(ref['clot_count'])].min(0);hi=ref['X'][:int(ref['clot_count'])].max(0);scale=250e-6
for i in range(-2,6):
 for j in range(-2,4):
  left=np.array([-375e-6+i*scale,-125e-6+j*scale]);right=left+scale
  a=np.maximum(left,lo[:2]);b=np.minimum(right,hi[:2])
  if np.all(b>a):cells.append(dict(i=i,j=j,center_m=(left+right)/2,nominal_size_m=[scale,scale],clipped_lower_m=a,clipped_upper_m=b,clipped_area_m2=float(np.prod(b-a)),rotation_sign=(-1)**(i+j),center_inside=bool(np.all((left+right)/2>=lo[:2]) and np.all((left+right)/2<=hi[:2]))))
dump(O/'SHEAR_PARAMETERIZATION.json',dict(modes=CASES,label=LABEL,counter_rotating=dict(center_m=[0,0,375e-6],width_m=125e-6,separation_m=300e-6,angle_rad=0,spatial_formula='Exactly prior double-Gaussian n cross grad psi, same C2 lab window'),checkerboard=dict(label='CHECKERBOARD_SHEAR_MECHANICAL_PROXY',scale_m=scale,actual_top_span_m=hi[:2]-lo[:2],cells_intersecting_top=len(cells),cell_centers_inside_top=sum(c['center_inside'] for c in cells),cells=cells,streamfunction='psi=sin(pi*(a+375um)/250um)*sin(pi*(b+125um)/250um)/(pi/250um)',normalization='Analytic |gradient psi| <= 1, continuous across cell boundaries. Edge cells honestly clipped; no specimen enlargement'),time=dict(rise_s=RISE,hold_s=HOLD,off_s=OFF,end_s=END),controller=SPEC,numerical=dict(BASE=BASE,MIN_H=MIN_H,max_dd=.1,energy_gate=.01),old_A='Work budget only, not same-surface-physics control'))
dump(O/'P0_SCOPE.json',dict(isolation=str(A),historical_evidence_read_only=True,deleted_files=[],reason='No confirmed unused disposable code identified; keep prior research intact and isolate all new paths.',git_sha=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip()))
(O/'provenance/git_existing.diff').write_bytes(subprocess.check_output(['git','-c',f'safe.directory={R}','diff'],cwd=R))
print('PREPARED',TASK,RUN,len(surf.tet),len(surf.tri),'checker cells',len(cells),flush=True)
