from common import *
from campaign_recovered import execute
source=O/'cases/surface_shear_1000pa_continuation_2/SUMMARY.json'
z=json.loads(source.read_text());assert z['status']=='COMPLETED'
peak=float(np.sqrt(100.*1000.))
dump(O/'CHUNK_MIDPOINT_DECISION.json',dict(peak_pa=peak,mode='SURFACE_SHEAR',basis_sha256=sha(source),reason='100 Pa completed without detachment; 1000 Pa formed a large cohesive component and smaller secondary debris. One intermediate case classifies the detachment/fragmentation regime. It is not bisection on maximum fragment size and assumes no size monotonicity.',additional_physical_cases=1))
execute('SURFACE_SHEAR',peak)
