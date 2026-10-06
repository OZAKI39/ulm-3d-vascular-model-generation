"""Re-execute only the two retained D2 cells with the original native harness."""
from common import *
import numpy as np
import subprocess

def main():
    old=read(PREVIOUS/'P2_CURVED_GEOMETRY_AUDIT.json')
    data=np.load(PREVIOUS/'P2_CURVED_GEOMETRY_SAMPLES.npz')
    exe=PREVIOUS/'runtime/native_geometry'
    assert sha(exe)==old['harness_binary_sha256']
    rt=HERE/'runtime';rt.mkdir(exist_ok=True)
    rows=[]
    for f in old['fixtures']:
        for variant in ['affine','curved']:
            tag=f['label']+'_'+variant
            np.savetxt(rt/(tag+'_nodes.txt'),data[tag+'_nodes'],fmt='%.17g')
            out=rt/(tag+'_baseline.json')
            subprocess.run([str(exe),str(rt/(tag+'_nodes.txt')),str(rt/(tag+'_baseline.txt')),str(out)],env=native_env(),check=True)
            actual=np.array(read(out)['native_construct_fluid_Rp']);prev=np.array(f['variants'][variant]['native_Rp'])
            common=data[tag+'_common_Rp']
            match=float(np.max(abs(actual-prev))/np.max(abs(prev)))
            mismatch=float(np.linalg.norm(actual-common)/np.linalg.norm(common))
            assert match<=1e-11
            assert (mismatch<=1e-11 if variant=='affine' else mismatch>1e-5)
            rows.append(dict(label=f['label'],cell_id=f['original_cell_id'],variant=variant,native_Rp=actual.tolist(),
                common_reference_Rp=common.tolist(),previous_result_relative_max_error=match,
                native_vs_common_relative_L2=mismatch))
    put('CURVED_MIXED_BLOCKER_BASELINE.json',dict(status='BLOCKER_REPRODUCED',fixtures=rows,
        baseline_binary_sha256=sha(exe),protocol_sha256=sha(HERE/'COMMON_GEOMETRY_HARNESS_PROTOCOL.json'),
        prior_A_to_D1_reexecuted=False,linear_pressure_extension='Will be recorded in native instrumented baseline before patched validation',
        CFD_steps=0,linear_solves=0))
    print('BLOCKER_REPRODUCED',[(r['label'],r['variant'],r['native_vs_common_relative_L2']) for r in rows])

if __name__=='__main__':main()
