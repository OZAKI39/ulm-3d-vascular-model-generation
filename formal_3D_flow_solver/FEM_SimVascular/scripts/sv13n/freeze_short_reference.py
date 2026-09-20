"""Freeze the already-complete step10; do not rerun the old solver."""
raise SystemExit('DEFERRED: latest user strategy retains old outputs only as development records')
import copy,json,subprocess,sys,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv13n import flow_gate
R=ROOT/'reports/sv1_3n';O=ROOT/'outputs/sv1_3n'
N=json.loads((ROOT/'configs/sv1_3n/convergence_strategy.json').read_text())['SHORT_REGRESSION_STEP']
original=json.loads((R/'OLD_PETSC_CPU_PROOF_20_acceptance.json').read_text());assert original['accepted']
d=copy.deepcopy(original);results=[f for f in d['results'] if Path(f['path']).stem==f'result_{N:03d}'];assert len(results)==1
path=O/Path(results[0]['path']).relative_to('outputs')
reload=json.loads(subprocess.check_output([sys.executable,'-B',ROOT/'scripts/sv13n/reload_field.py',path],text=True))
source=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
m=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',source['accepted_solution']['Q_target_m3_s'],source['policy']['Umean_m_s'])
u,p=m.read(path);measured=m.measure(u,p)
rows=[r for r in d['history']['linear_solves'] if r['step']<=N];assert len({r['step'] for r in rows})==N
d.update(name='OLD_CPU_SHORT_REFERENCE',origin_case=original['name'],results=results,reload=[reload],VTU_count=1,
 steps_completed=N,checkpoint_step=N,requested_steps=N,measurement=measured,linear_solves=len(rows),runtime_semantics=d['runtime_semantics'][:len(rows)],
 wall_time_s=None,observed_checkpoint_elapsed_s=rows[-1]['reported_elapsed_s'],runtime_scope='Solver-reported elapsed at final converged step10 row; rounded log value, not isolated process wall time.',
 original_safe_stopped_wall_time_s=original['wall_time_s'],original_completed_steps=original['steps_completed'],steady_reference=False)
d['history']['linear_solves']=rows
d['velocity_finite']=reload['velocity_finite'];d['pressure_finite']=reload['pressure_finite'];d['wall_noslip_pass']=measured['wall_noslip_pass'];d['flows_finite']=all(math.isfinite(v) for v in [measured['Q_in_m3_s'],*measured['outlet_flows_m3_s'].values(),measured['epsilon_mass']])
flow_gate(d,gpu=False,steps=N)
(R/'OLD_CPU_SHORT_REFERENCE_acceptance.json').write_text(json.dumps(d,indent=2)+'\n')
print('OLD_CPU_SHORT_REFERENCE: matched step10 finite fields/flows/wall/reload PASS; not a steady reference.')
