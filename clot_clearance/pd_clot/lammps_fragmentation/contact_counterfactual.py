"""Single-case contact interventions; original map files and physics stay frozen."""
from pathlib import Path
import json, shutil, difflib, time
import numpy as np
from .config import ROOT, dump
from . import runner
from .failure_map import analyze_case

O=ROOT/'verification/lammps_contact_cascade_diagnosis_and_minimal_fix_v1'
R=ROOT/'results/lammps_contact_cascade_diagnosis_and_minimal_fix_v1'
OLD=ROOT/'verification/lammps_coherent_failure_force_time_map_v1'
BASE=OLD/'cases/11_LOCAL_GAUSSIAN_SHEAR_A100_D2'

def geometry():
    z=np.load(OLD/'INTERFACE_GEOMETRY.npz');g={k:z[k] for k in z.files}
    g['clot_count']=int(g['clot_count']);g['info']=json.loads((OLD/'INTERFACE_GEOMETRY.json').read_text())
    return g

def intervene(text, model):
    if model=='PAIR_ZERO_DIAGNOSTIC':
        source='pair_style bpm/peri\npair_coeff * * 3.7226771877469815e+20 0.0001505\n'
        target='pair_style zero 0.0001505\npair_coeff * *\n'
    elif model=='MODIFIED_BPM_PERI':
        source='pair_style bpm/peri\n';target='pair_style bpm/peri contact_distance_factor 1.0\n'
    else:raise ValueError(model)
    assert text.count(source)==1
    changed=text.replace(source,target);assert changed.replace(target,source)==text
    return changed

def prepare(model):
    name={'PAIR_ZERO_DIAGNOSTIC':'CONTACT_OFF','MODIFIED_BPM_PERI':'CONTACT_MODIFIED'}[model]
    p=O/name;p.mkdir(exist_ok=False)
    for f in ['body.data','LOAD_PROVIDER.json']:shutil.copy2(BASE/f,p/f)
    original=(BASE/'in.lammps').read_text();changed=intervene(original,model);(p/'in.lammps').write_text(changed)
    diff=''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile='CONTACT_ON/in.lammps',tofile=name+'/in.lammps'))
    (p/'input.diff').write_text(diff)
    q=json.loads((BASE/'request.json').read_text());q.update(case_id=name,contact_model=model,baseline_case_id=BASE.name,contact_scope='DIAGNOSTIC_ONLY' if model=='PAIR_ZERO_DIAGNOSTIC' else 'SINGLE_ENGINEERING_CANDIDATE')
    dump(p/'request.json',q)
    with (O/'INPUT_DIFF.md').open('a') as f:f.write('\n## '+name+'\n\nOnly the displayed pair-contact command changes. Initial data are byte-identical; every other input line, integration/fix order, communication cutoff, neighbor skin, MPI ranks and saved step is retained.\n\n```diff\n'+diff+'```\n')
    return dict(path=p,request=q,geometry=geometry(),load=json.loads((p/'LOAD_PROVIDER.json').read_text()),steps=q['explicit_steps'])

def run(model, executable=None):
    case=prepare(model);start=time.perf_counter();saved=runner.LMP
    try:
        if executable is not None:runner.LMP=Path(executable)
        runner.run_lammps(case)
    finally:runner.LMP=saved
    result=analyze_case(case);result.update(contact_model=model,contact_scope=case['request']['contact_scope'],baseline_case_id=BASE.name,run_and_postprocess_wall_s=time.perf_counter()-start)
    dump(case['path']/'result.json',result)
    print('CONTACT_RESULT',model,result['failure_mode'],result['bulk_broken_fraction'],result['debris_fraction'],flush=True)
    return result
