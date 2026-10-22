"""Freeze all inputs and unseen seeds before F10AE physics."""
import h5py, numpy as np, subprocess, platform, sys
from scripts.straight_tube_exp1.f10ab.prepare import shape
from scripts.straight_tube_exp1.f10ab.metrics import geometry
from scripts.straight_tube_exp1.f10s.common import context
from .common import *
def main():
    BASE.mkdir(parents=True,exist_ok=True)
    assert not (BASE/'PREREGISTRATION.json').exists()
    cfg=read(CONFIG)
    with h5py.File(AC/'INITIAL_STATES.h5') as old,h5py.File(BASE/'INITIAL_STATES.h5','x') as h:
        h['S2_r0']=old['S2_r0'][:]
        for label,seeds in [('validation',cfg['validation_seeds']),('test',cfg['final_test_seeds'])]:
            for i,seed in enumerate(seeds):h[f'S2_{label}_{i}']=shape('S2',seed)
        audit={k:geometry(h[k][:]) for k in h}
        c=context()
        for k in h:
            x=h[k][:]
            assert len(x)==64 and audit[k]['proximity_components']==1
            assert min(c['cavity'].closest_wall(q)[0]-R for q in x)>100e-6
    write(BASE/'INITIAL_STATE_AUDIT.json',audit)
    protected={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'src/mevo_mb').rglob('*.py')}
    for p in [CONFIG,BASE/'INITIAL_STATES.h5',OLD/'inputs/response.h5',OLD/'cache/coarse.h5',OLD/'cache/fine.h5',AC/'INITIAL_STATES.h5',AC/'FINAL_GATE.json',AD/'2400000/response_shell_1.0.h5',AD/'2600000/response_shell_1.0.h5']:
        protected[str(p.relative_to(ROOT))]=sha(p)
    # All raw historical trajectories protected, including failures.
    history={str(p.resolve()):sha(p) for parent in [AC,AB] for p in (parent/'runs').rglob('trajectory.h5')}
    write(BASE/'HISTORICAL_RAW_SHA256.json',history)
    write(BASE/'PREREGISTRATION.json',dict(identity=identity(),config=cfg,protected_sha256=protected,historical_raw_sha256=sha(BASE/'HISTORICAL_RAW_SHA256.json'),independent_test_policy='Freeze selected actions before loading test trajectories; no action adjustment after final tests'))
    write(BASE/'ENVIRONMENT.json',dict(python=sys.version,platform=platform.platform(),packages=subprocess.check_output(['uv','pip','freeze','--python',sys.executable],text=True),nvidia=subprocess.check_output(['nvidia-smi'],text=True)))
    print('F10AE preregistered',len(protected),len(history),flush=True)
if __name__=='__main__':main()
