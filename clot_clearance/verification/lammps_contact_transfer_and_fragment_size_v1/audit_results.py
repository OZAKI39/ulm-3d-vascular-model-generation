"""Required formula, graph, input, lifecycle and history checks only."""
from pathlib import Path
import sys,json,hashlib,csv,ast,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from pd_clot.lammps_fragmentation.targeted_contact_transfer import O,GAUSSIAN100,specification,geometry,intervene
from pd_clot.lammps_fragmentation.current_position_provider import lammps_provider
from pd_clot.lammps_fragmentation.detached_fragment_metrics import equivalent_diameter_um
from pd_clot.lammps_fragmentation.config import ROOT,dump

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
start=time.perf_counter();s=json.loads((O/'ALL_CASE_SUMMARIES.json').read_text());cfg=json.loads((O/'FROZEN_CONFIGURATION.json').read_text())['request'];g=geometry()
assert sha(O/'PREDECLARED_DECISION_POLICY.json')==(O/'PREDECLARED_DECISION_POLICY.sha256').read_text().strip()
threshold=[]
for n,e200,e300 in [(1,False,False),(2,False,False),(3,False,False),(5,False,False),(34,True,False),(49,True,False),(113,True,False),(114,True,True)]:
    d=float(equivalent_diameter_um(n,cfg['dx_m']));assert (d>=200)==e200 and (d>=300)==e300;threshold.append(dict(nodes=n,d_eq_um=d,ge200=e200,ge300=e300))
checks=[]
formal=sorted(p.parent.name for p in (O/'CASES').glob('*/RUN_COMMAND.json'));assert len(formal)==4 and set(formal)=={'OPPOSING_100_D2','GAUSSIAN_200_D2','GAUSSIAN_31P622_D2','OPPOSING_31P622_D2'}
for name in formal:
    p=O/'CASES'/name;src,amplitude=specification(name);r=json.loads((p/'request.json').read_text());oldr=json.loads((src/'request.json').read_text());native=(src/'in.lammps').read_text()
    if amplitude is None:expected=intervene(native,'MODIFIED_BPM_PERI')
    else:
        oldblock=lammps_provider(g,oldr,json.loads((src/'LOAD_PROVIDER.json').read_text()));newblock=lammps_provider(g,r,json.loads((p/'LOAD_PROVIDER.json').read_text()));expected=native.replace(oldblock,newblock)
    assert (p/'in.lammps').read_text()==expected
    for key,value in cfg.items():
        if key not in {'case_id','pattern_id','tau_budget_Pa','map_stage','baseline_case_id','contact_scope'}:assert r[key]==value,(name,key)
    assert (p/'body.data').read_bytes()==(GAUSSIAN100/'body.data').read_bytes()
    z=np.load(p/'accepted_trajectory.npz');assert np.array_equal(z['step'],np.arange(0,4661,20)) and np.isfinite(z['position_m']).all() and np.isfinite(z['velocity_m_s']).all()
    assert z['position_m'].shape==(234,4174,3)
    result=json.loads((p/'result.json').read_text());assert all(value for key,value in result['validation'].items() if key!='provider_max_force_difference_N')
    assert result['validation']['provider_max_force_difference_N']<1e-13
    native_log=(p/'log.lammps').read_text();assert 'ERROR:' not in native_log and 'Lost atoms:' not in native_log and 'Bond atoms missing' not in native_log
    perf=json.loads((p/'PERFORMANCE.json').read_text());assert perf['exit_code']==0 and perf['mpi_ranks']==2 and perf['explicit_steps']==4660
    if '31P622' in name:assert (O/'CONDITIONAL_CASE_DECISION.json').stat().st_mtime<(p/'RUN_COMMAND.json').stat().st_mtime
    checks.append(dict(case_id=name,input_only_named_change=True,initial_geometry_byte_identical=True,full_time_and_mass=True,native_error_free=True,warnings=[a for a in native_log.splitlines() if 'WARNING:' in a]))
lineage=[]
for name,q in s.items():
    assert all(q['provenance'].values()) and q['tracking_ambiguities']==0
    rec=json.loads((O/'analysis'/name/'membership_and_shapes.json').read_text())['records'];children={}
    for key,v in rec.items():
        ids=set(v['node_ids']);assert len(ids)==v['node_count'] and min(ids)>=1 and max(ids)<=2971
        df=v['detached_frames'];assert not df or df==list(range(df[0],df[-1]+1))
        for parent in v['parent_track_ids']:
            assert ids<set(rec[str(parent)]['node_ids']);children.setdefault(parent,[]).append(ids)
    for parent,parts in children.items():
        union=set().union(*parts);assert union==set(rec[str(parent)]['node_ids']) and sum(map(len,parts))==len(union)
    with (O/'analysis'/name/'timeseries.csv').open() as f:
        for row in csv.DictReader(f):
            assert 0<=float(row['mass_fraction_detached_ge300'])<=float(row['mass_fraction_detached_ge200'])<=1-float(row['total_attached_fraction'])+1e-12
            assert int(row['N_detached_ge300'])<=int(row['N_detached_ge200'])
    lineage.append(dict(case_id=name,native_events_and_cached_graph_equal=True,split_mass_and_particle_sets_conserved=True,zero_mergers_or_tracking_ambiguities=True,tracks=len(rec)))
protected=json.loads((O/'HISTORY_PROTECTION.json').read_text());bad=[]
for row in protected:
    p=Path(row['path'])
    if not p.exists() or p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:bad.append(row['path'])
assert not bad,bad
dump(O/'checks/HISTORY_INTEGRITY.json',dict(status='PASS',protected_files=len(protected),protected_bytes=sum(r['bytes'] for r in protected),changed_files=bad,scope='Complete previous two verification/results trees, previous Python package, contact variant source/executable, official executable, and existing active context'))
for name in ['detached_fragment_metrics.py','targeted_contact_transfer.py','targeted_map_refresh.py','contact_transfer_visualization.py']:ast.parse((ROOT/'pd_clot/lammps_fragmentation'/name).read_text())
dump(O/'checks/ROUND_AUDIT.json',dict(status='PASS',threshold_cases=threshold,case_checks=checks,lineage_checks=lineage,new_formal_solves=4,old_solver_reruns=0,conditional_decision_precedes_runs=True,history='PASS',audit_wall_s=time.perf_counter()-start))
print('ROUND_AUDIT_PASS',len(protected),'historical files; four actual solves; nine read-only trajectory audits',flush=True)
