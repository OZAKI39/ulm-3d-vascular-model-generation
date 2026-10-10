"""Bounded transfer campaign, reusing the verified contact variant unchanged."""
from pathlib import Path
import json,shutil,difflib,time
import numpy as np
from .config import ROOT,dump
from . import runner
from .contact_counterfactual import geometry,intervene,OLD,BASE
from .current_position_provider import prepare_provider,lammps_provider
from .failure_map import analyze_case,csv_out
from .detached_fragment_metrics import analyze,focused_checks

O=ROOT/'verification/lammps_contact_transfer_and_fragment_size_v1'
R=ROOT/'results/lammps_contact_transfer_and_fragment_size_v1'
PREVIOUS=ROOT/'verification/lammps_contact_cascade_diagnosis_and_minimal_fix_v1'
GAUSSIAN100=PREVIOUS/'CONTACT_MODIFIED'
OPPOSING_NATIVE=OLD/'cases/14_OPPOSING_SHEAR_A100_D2'
EXE=ROOT/'tools/lammps-contact-onset-v1-build/lmp'
REFERENCES={'GAUSSIAN_100_MODIFIED':(GAUSSIAN100,1.0),'GAUSSIAN_100_NATIVE':(BASE,1.35),'OPPOSING_100_NATIVE':(OPPOSING_NATIVE,1.35)}

def stage_zero():
    assert (O/'HISTORY_PROTECTION.json').exists()
    dump(O/'checks/SIZE_AND_ATTACHMENT_TESTS.json',focused_checks())
    result=analyze(GAUSSIAN100,O/'analysis/GAUSSIAN_100_MODIFIED','GAUSSIAN_100_MODIFIED',geometry(),1.0)
    membership=json.loads((O/'analysis/GAUSSIAN_100_MODIFIED/membership_and_shapes.json').read_text())['records']
    hits=[v for v in membership.values() if v['node_count']==49 and v['detached_frames'] and v['detached_frames'][0]==50]
    assert len(hits)==1,hits
    life=np.genfromtxt(O/'analysis/GAUSSIAN_100_MODIFIED/lifetimes.csv',delimiter=',',names=True,dtype=None,encoding='utf8')
    z=life[life['track_id']==hits[0]['track_id']][0]
    example={k:z[k].item() for k in life.dtype.names};example['node_ids']=hits[0]['node_ids'];example['crosscheck']='49-node DETACHED component independently recovered at saved step 1000 / 0.500 ms; all native events and interface bonds reread'
    dump(O/'HISTORICAL_GAUSSIAN_REANALYSIS.json',dict(summary=result,example_49_nodes=example))
    for name,(path,factor) in REFERENCES.items():
        if name!='GAUSSIAN_100_MODIFIED':analyze(path,O/'analysis'/name,name,geometry(),factor)
    print('STAGE_0_PASS',example,flush=True)

def specification(name):
    if name=='OPPOSING_100_D2':return OPPOSING_NATIVE,None
    if name=='GAUSSIAN_200_D2':return GAUSSIAN100,200.
    if name=='GAUSSIAN_31P622_D2':return OLD/'cases/31_LOCAL_GAUSSIAN_SHEAR_A31.62278_D2',None
    if name=='OPPOSING_31P622_D2':return OLD/'cases/32_OPPOSING_SHEAR_A31.62278_D2',None
    raise ValueError(name)

def prepare(name):
    src,amp=specification(name);r=json.loads((src/'request.json').read_text());oldr=r.copy();g=geometry();original=(src/'in.lammps').read_text();oldload=json.loads((src/'LOAD_PROVIDER.json').read_text())
    if amp is None:
        changed=intervene(original,'MODIFIED_BPM_PERI');load=oldload
    else:
        r['tau_budget_Pa']=amp;load=prepare_provider(g,r);oldblock=lammps_provider(g,oldr,oldload);newblock=lammps_provider(g,r,load);assert original.count(oldblock)==1
        changed=original.replace(oldblock,newblock)
        diff=[x for x in difflib.ndiff(original.splitlines(),changed.splitlines()) if x.startswith(('- ','+ '))]
        assert len(diff)==6 and all(x[2:].startswith('variable load_f') for x in diff)
        assert changed.replace(newblock,oldblock)==original
    r.update(case_id=name,contact_model='MODIFIED_BPM_PERI',contact_model_id='MODIFIED_1P0DX',contact_distance_factor=1.,map_stage='TARGETED_CONTACT_TRANSFER',baseline_case_id=src.name,contact_scope='CONDITIONAL_EXPLORATORY_CANDIDATE',round_name='LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1')
    frozen=json.loads((O/'FROZEN_CONFIGURATION.json').read_text())['request'];allowed={'case_id','pattern_id','tau_budget_Pa','map_stage','baseline_case_id','contact_scope'}
    for key,value in frozen.items():
        if key not in allowed:assert r[key]==value,(name,key,r[key],value)
    path=O/'CASES'/name;path.mkdir(exist_ok=False);shutil.copy2(src/'body.data',path/'body.data');assert (path/'body.data').read_bytes()==(GAUSSIAN100/'body.data').read_bytes()
    dump(path/'request.json',r);dump(path/'LOAD_PROVIDER.json',load);(path/'in.lammps').write_text(changed)
    delta=''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile=str(src/'in.lammps'),tofile=name+'/in.lammps'));(path/'input.diff').write_text(delta)
    dump(path/'INPUT_PROVENANCE.json',dict(source_path=str(src),allowed_physical_change='tau_budget_Pa: 100 -> 200 only' if amp else 'native pair onset 1.35 -> explicit 1.0 only',initial_data_byte_identical=True,all_other_physical_request_fields_equal=True,nominal_amplitude_Pa=r['tau_budget_Pa'],actual_pulse_duration_s=r['duration_s'],dt_s=r['dt_s'],steps=r['explicit_steps'],dump_every=r['dump_every'],mpi_ranks=r['mpi_ranks']))
    return dict(path=path,request=r,geometry=g,load=load,steps=r['explicit_steps'])

def run_one(name):
    if '31P622' in name:
        decision=json.loads((O/'CONDITIONAL_CASE_DECISION.json').read_text());assert name in decision['selected'] and decision['allow']
    assert len(list((O/'CASES').glob('*/RUN_COMMAND.json')))<4
    case=prepare(name);oldexe=runner.LMP
    try:runner.LMP=EXE;runner.run_lammps(case)
    finally:runner.LMP=oldexe
    result=analyze_case(case);result.update(contact_model_id='MODIFIED_1P0DX',contact_distance_factor=1.0,scope='EXPLORATORY_ONLY');dump(case['path']/'result.json',result)
    summary=analyze(case['path'],O/'analysis'/name,name,case['geometry'],1.)
    dump(case['path']/'FRAGMENT_SIZE_RESULT.json',summary)
    print('FORMAL_CASE_COMPLETE',name,summary['final'],flush=True);return summary

def stage_one_two():
    assert (O/'HISTORICAL_GAUSSIAN_REANALYSIS.json').exists()
    for name in ['OPPOSING_100_D2','GAUSSIAN_200_D2']:run_one(name)

def decide_conditional():
    policy=json.loads((O/'PREDECLARED_DECISION_POLICY.json').read_text())
    def read(name):return json.loads((O/'analysis'/name/'summary.json').read_text())
    new=read('OPPOSING_100_D2');old=read('OPPOSING_100_NATIVE');strong=read('GAUSSIAN_200_D2')
    def normal(q):return q['n_steps']==4660 and q['time_end_s']==.00233 and not q['sustained_severe_overlap'] and q['tracking_ambiguities']==0 and all(q['provenance'].values())
    a=np.load(OPPOSING_NATIVE/'accepted_trajectory.npz');b=np.load(O/'CASES/OPPOSING_100_D2/accepted_trajectory.npz');t=a['step']*.5e-6
    first=min(x for q in [new,old] for x in [q['first_bulk_break_time_s'],q['first_interface_break_time_s']] if x is not None);pre=t<first
    pe=float(np.linalg.norm(a['position_m'][pre]-b['position_m'][pre],axis=2).max()/50e-6);ve=float(np.linalg.norm(a['velocity_m_s'][pre]-b['velocity_m_s'][pre],axis=2).max()/np.sqrt(.94))
    parity=pe<=policy['prebreak_position_tolerance_dx'] and ve<=policy['prebreak_velocity_tolerance_cS'];drop=old['final']['debris_fraction']-new['final']['debris_fraction']
    supported=normal(new) and parity and new['max_active_contact_pairs']>0 and new['final']['failure_mode']!='PULVERIZATION' and drop>=policy['significant_terminal_debris_drop']
    decision=dict(status='SUPPORTED_FOR_TESTED_PATTERNS' if supported else 'NOT_ESTABLISHED',normal=normal(new),prebreak_parity=parity,max_prebreak_position_difference_dx=pe,max_prebreak_velocity_difference_cS=ve,prebreak_saved_frames=int(pre.sum()),terminal_debris_absolute_drop=drop,new_terminal_mode=new['final']['failure_mode'],native_terminal_mode=old['final']['failure_mode'],contact_active_peak=new['max_active_contact_pairs'],sustained_severe_overlap=new['sustained_severe_overlap'],scope='Only tested Gaussian and Opposing at 100 Pa / 2Tcross; 200 Pa response assessed separately')
    dump(O/'CROSS_PATTERN_CONTACT_DECISION.json',decision)
    allow=supported and normal(strong)
    selected=['GAUSSIAN_31P622_D2','OPPOSING_31P622_D2'] if allow else []
    frozen=[]
    for name in selected:
        p,_=specification(name);request=json.loads((p/'request.json').read_text());frozen.append(dict(case_id=name,source=str(p),amplitude_Pa=request['tau_budget_Pa'],pattern=request['pattern_id'],contact_factor=1.,pulse_duration_s=request['duration_s'],dt_s=request['dt_s'],end_step=request['explicit_steps'],sole_input_change='pair_style bpm/peri contact_distance_factor 1.0'))
    conditional=dict(allow=allow,selected=selected,rejected=[] if allow else ['GAUSSIAN_31P622_D2','OPPOSING_31P622_D2'],mandatory_cases_normal=normal(new) and normal(strong),cross_transfer_supported=supported,reason='New 100 Pa modes differ from native pulverization; the two existing 31.6228 Pa native points bracket different pattern transitions, so same-amplitude modified points are informative without interpolating a full map.' if allow else 'Mandatory normality/transfer gate not met; no additional mechanics authorized by this branch.',frozen_inputs=frozen,max_conditional_solves=2,written_before_conditional_solves=True)
    dump(O/'CONDITIONAL_CASE_DECISION.json',conditional);print('CONDITIONAL_DECISION',conditional,flush=True);return conditional

def stage_three():
    decision=json.loads((O/'CONDITIONAL_CASE_DECISION.json').read_text())
    for name in decision['selected']:
        run_one(name);src,_=specification(name);analyze(src,O/'analysis'/(name+'_NATIVE'),name+'_NATIVE',geometry(),1.35)
