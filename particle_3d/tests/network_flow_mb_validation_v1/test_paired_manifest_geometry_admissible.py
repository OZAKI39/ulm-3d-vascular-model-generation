import numpy as np
from particle_3d.injection_admission import FiniteSizeAdmission
def test_paired_manifest_geometry_admissible(events,oldenv,newenv):
    for e in events:
        checks=[]
        for env in [oldenv,newenv]:
            p,s,d=FiniteSizeAdmission(None,None,wall=env.wall,field=env.field).check(e,np.array(e['birth_center_m']),{})
            checks.append((p is not None,s,d,env.field.locate(e['birth_center_m'])[0]))
        assert checks[0]==checks[1] and checks[0][0] and checks[0][3]>=0
