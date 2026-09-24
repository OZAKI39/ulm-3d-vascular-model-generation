import numpy as np
from particle_3d.rbc_orientation import short_axis
def test_haar(source):
 q=source.orientation(100000); p=np.array([short_axis(x) for x in q])
 assert abs(np.linalg.norm(q,axis=1)-1).max()<1e-14
 assert abs(p.mean(0)).max()<.007
 assert abs((p*p).mean(0)-1/3).max()<.006
