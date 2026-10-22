"""Freeze source, historical identities, bounded hypotheses and legal source order."""
import time, json, subprocess, importlib.metadata as md
import numpy as np
from .common import *
def main():
    start=time.perf_counter();BASE.mkdir(parents=True,exist_ok=True)
    assert not (BASE/'PREREGISTRATION.json').exists()
    history=read(AE/'INPUT_OUTPUT_SHA256.json');g=read(AE/'FINAL_GATE.json')['gates']
    assert g['F10AA_BACKEND_VERIFIED']['status']=='PASS' and g['STEP_SIZE_CONVERGENCE_PASS']['status']=='PASS'
    protected={str(p.relative_to(ROOT)):sha(p) for folder in ['src/mevo_mb/fast_sim','src/mevo_mb/acoustics','src/mevo_mb/control','src/mevo_mb/coupling','scripts/straight_tube_exp1/f10aa'] for p in (ROOT/folder).glob('*.py')}
    for p in [AE/'low_pressure/diagnostic_response.h5',OLD/'inputs/response.h5',OLD/'CONFIG_SNAPSHOT.json',AE/'INITIAL_STATES.h5',AE/'LOW_PRESSURE_ADMISSION.json']:
        protected[str(p.relative_to(ROOT))]=sha(p)
    for f in [2400000,2500000,2600000]:
        p=AD/str(f)/'response_shell_1.0.h5';protected[str(p.relative_to(ROOT))]=sha(p)
    write(BASE/'HISTORICAL_AE_SHA256.json',history)
    write(BASE/'PREREGISTRATION.json',dict(protocol=read(CONFIG),protected_sha256=protected,baseline_commit='306d899fa4a9be18b8a926688a6e25d389193138',historical_scientific_execution=read(AE/'WORKFLOW_COMPLETE.json')['identity'],identity=identity()))
    elements=read(OLD/'CONFIG_SNAPSHOT.json')['geometry']['pzt'];centers=np.array([e['center_local_m'] for e in elements]);order=np.argsort(np.mod(np.arctan2(centers[:,1],centers[:,0]),2*np.pi));phase=np.zeros(4);phase[order]=np.arange(4)*np.pi/2
    actions={}
    def add(name,a,p):actions[name]=dict(name=name,amplitudes=a,phases=list(p),frequency=2500000.,scope='Pa-equivalent sources, not hardware voltage')
    add('OFF',[0]*4,[0]*4);add('X_PAIR',[0,0,1000,1000],[0,0,0,0]);add('X_180',[0,0,1000,1000],[0,0,0,np.pi]);add('Y_PAIR',[1000,1000,0,0],[0]*4);add('XY_INPHASE',[1000]*4,[0]*4);add('XY_QUADRATURE',[1000]*4,[np.pi/2,np.pi/2,0,0]);add('ROT_PLUS',[1000]*4,phase);add('ROT_MINUS',[1000]*4,-phase);add('X_SINGLE',[0,0,1000,0],[0]*4)
    write(BASE/'ACTIONS.json',actions);write(BASE/'SOURCE_ORDER.json',dict(elements=elements,ccw_indices=order,ccw_names=[elements[i]['id'] for i in order],positive_phase=phase,negative_phase=-phase,equal_squared_amplitude_budget=True,scope='Static coherent phase winding, not a physically rotating amplitude pattern',identity=identity()))
    versions={n:md.version(n) for n in ['numpy','scipy','warp-lang','torch','botorch','ax-platform','smac','ConfigSpace']}
    resources=dict(gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu,driver_version','--format=csv'],text=True),versions=versions,wall_s=time.perf_counter()-start,identity=identity())
    write(BASE/'ENVIRONMENT.json',resources)
    print('PREPARED',versions,flush=True)
if __name__=='__main__':main()
