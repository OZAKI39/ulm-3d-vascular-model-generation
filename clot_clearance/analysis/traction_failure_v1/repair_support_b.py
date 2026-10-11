"""One explicitly frozen contact-rim graph repair before any force map exists."""
from analysis.traction_failure_v1.common import *
import shutil
assert not (B/'CONTACT_RIM_SUPPORT_REPAIR_FREEZE.json').exists()
old=json.loads((B/'COVERAGE.json').read_text());s=np.load(PREV/'refined/CLOT_SURFACE_TRACTIONS.npz');faces=old['unsupported_faces'];assert len(faces)==10
assert all(s['contact_rim_distance_m'][faces]<DX/4)
history=B/'initial_support_rejection';history.mkdir()
for name in ['FORMAL_MAPPING_START.json','COVERAGE.json','FINAL_STATUS.json']:
 shutil.move(B/name,history/name)
save(B/'CONTACT_RIM_SUPPORT_REPAIR_FREEZE.json',dict(time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),faces=faces,area_fraction=float(s['triangle_area_m2'][faces].sum()/s['triangle_area_m2'].sum()),reason='A tangent-plane halfspace is not a global interior test at a curved thin contact rim. All ten rejected faces already have local material support. Replace only their halfspace test by connected-surface face-center graph distance <= original150.5um horizon to the closest exposed triangle of each candidate node.',unchanged='Face-node Euclidean radius<=delta, material-to-surface distance<=dx, positive normal alignment, nonnegative weights; all original conservation/correction/direction/concentration gates.',algorithm='Actual triangle adjacency via shared native mesh edge; shortest path along centroid adjacency. Centroid chord-graph metric is not exact intrinsic geodesic or a guaranteed upper bound; links never bridge absent facets. Only this single repair, no radius/kernel/threshold scan.',previous_frozen_design_sha256=sha(B/'MAPPING_DESIGN_FREEZE.json'),prior_attempt='Stopped before initial force mapping; no integrated or candidate node forces existed',support_repair_attempts=1))
from analysis.traction_failure_v1.map_b import run
run()
