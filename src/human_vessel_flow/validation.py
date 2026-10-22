"""Read-only checks with explicit scope; copied data never acquires a new identity."""
from pathlib import Path
from collections import Counter
import ast, gzip, hashlib, json, sys, traceback
import numpy as np
import pyvista as pv
from .paths import ROOT, read_json, write_json, sha256, fresh_directory
from . import geometry, meshing, flow_field, wss, particles

def preservation(root=ROOT):
    manifest=read_json(root/'docs/provenance/SOURCE_MANIFEST_BEFORE_COPY.json')
    moves={r['original_destination']:r['destination'] for r in read_json(root/'docs/provenance/LAYOUT_MOVES.json')}
    checked=[]
    for row in manifest:
        if row['category'] in {'source','test'}:continue
        relative=moves.get(row['destination'],row['destination'])
        p=root/relative
        assert p.is_file() and not p.is_symlink(), relative
        assert sha256(p)==row['sha256'],relative
        checked.append(relative)
    return dict(status='PASS',immutable_files=len(checked),checked=checked)

def cfd(root=ROOT):
    result={}
    for name in ['balanced_pressure_final','calibration_equal_split']:
        field=flow_field.load(name,root);metric=flow_field.measurements(field,root)
        assert metric['epsilon_mass']<=field.policy['mass_limit'] and metric['wall_noslip_pass']
        assert abs(metric['Q_in_m3_s']-3e-7)/3e-7<1e-12
        surface,raw=wss.calculate(field)
        frozen=pv.read(root/'cases'/name/'postprocess/wall_wss_si.vtp')
        assert np.array_equal(surface['Parent_tetra_zero_based'],frozen['Parent_tetra_zero_based'])
        np.testing.assert_allclose(surface['WSS_raw_Pa'],frozen['WSS_raw_Pa'],rtol=1e-12,atol=1e-12)
        expected=read_json(root/'cases'/name/'postprocess/INDEPENDENT_CHECK.json')
        np.testing.assert_allclose([metric['outlet_flows_m3_s'][k] for k in sorted(metric['outlet_flows_m3_s'])],
            [expected['measurements']['outlet_flows_m3_s'][k] for k in sorted(metric['outlet_flows_m3_s'])],rtol=1e-12,atol=0)
        # Final native checkpoint is a historical record, not a new restart test.
        from flow_solver_support.transient_checks import checkpoint_audit
        case=root/'cases'/name;execution=read_json(case/'reports/execution.json')
        cp=case/execution['final_vtu'];cp=cp.with_name('stFile_%03d.bin'%execution['final_step'])
        checkpoint=checkpoint_audit(cp,execution['final_step'],field.policy['dt_s'],ranks=execution['MPI_ranks'])
        result[name]=dict(status='PASS_READ_ONLY',flow_sha256=sha256(root/field.record['flow']),
            points=field.mesh.n_points,tetra=field.mesh.n_cells,correspondence=field.correspondence,
            measurements=metric,raw_WSS_faces=len(raw['magnitude']),exact_zero_WSS_faces=int((raw['magnitude']==0).sum()),
            raw_WSS_area_mean_Pa=float(np.average(raw['magnitude'],weights=raw['area'])),checkpoint=checkpoint,
            role=field.record['physical_role'],local_physics_certification=False)
    return dict(status='PASS',cases=result)

def population(root=ROOT):
    run=root/'results/population/realization_0000_contactfix'
    manifest=read_json(run/'MANIFEST.json');summary=particles.historical_summary(root)
    sources=list(particles.ledger(run/'source.jsonl.gz'));terminals=list(particles.ledger(run/'terminals.jsonl.gz'))
    assert len(sources)==len(terminals)==1565
    assert [s['particle_id'] for s in sources]==list(range(1,1566))
    assert sorted(t['particle_id'] for t in terminals)==list(range(1,1566))
    assert all(s['status']=='ACCEPTED_BIRTH' for s in sources)
    counts=Counter(t['status'] for t in terminals)
    outlets=Counter(t['outlet'] for t in terminals if t['status']=='EXITED')
    assert dict(outlets)=={'OUTLET_01':458,'OUTLET_02':527,'OUTLET_03':437}
    assert counts['SOLVER_FAILURE']==2 and sum(counts.values())==1565
    assert summary['censored']==141 and summary['scientific_integrity_pass'] is False
    births={s['particle_id']:s for s in sources};last=-np.inf;rows=0;times=0;chunks=[];observed=set()
    for record in manifest['chunks']:
        p=run/record['path'];assert sha256(p)==record['file_sha256']
        with np.load(p) as a:
            t=a['state_time_s'];off=a['state_offsets'];ids=a['particle_id']
            assert len(ids)==record['row_count'] and len(t)==record['state_count']
            assert t[0]>=last and np.all(np.diff(t)>=0)
            assert off[0]==0 and off[-1]==len(ids) and np.all(np.diff(off)>=0)
            repeated=np.repeat(t,np.diff(off));bt=np.array([births[int(i)]['birth_time_s'] for i in ids])
            assert np.array_equal(a['birth_time'],bt)
            assert np.array_equal(a['source_event_id'],ids) and np.all(a['realization_id']==0)
            np.testing.assert_allclose(a['age'],repeated-bt,rtol=0,atol=2e-15)
            for k in ['x','v','omega','q','radius','diameter','birth_time','age']:
                assert np.isfinite(a[k]).all(),(p.name,k)
            assert np.all(a['radius']>0) and np.all(a['age']>=0)
            observed.update(map(int,np.unique(ids)));rows+=len(ids);times+=len(t);last=float(t[-1])
            chunks.append(dict(path=record['path'],sha256=record['file_sha256'],rows=len(ids),states=len(t)))
    assert observed==set(births) and rows==manifest['recorded_particle_states']
    ledger_hashes={}
    for name,expected in manifest['ledger_scientific_hashes'].items():
        digest=hashlib.sha256()
        with gzip.open(run/(name+'.jsonl.gz'),'rb') as f:
            for block in iter(lambda:f.read(2**20),b''):digest.update(block)
        assert digest.hexdigest()==expected,name
        ledger_hashes[name]=expected
    return dict(status='PASS_PRESERVED_ACCOUNTING',summary=summary,terminal_counts=dict(counts),
        particle_rows=rows,snapshot_count=times,chunks=chunks,ledger_scientific_hashes=ledger_hashes,
        all_trajectories_successful=False,new_fixed_window_formal_run=False)

def scientific_gates(root=ROOT):
    ref=root/'diagnostics/reference_pipe'
    expected={'FIELD_VALIDATION.json':'REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC',
              'V2_FIELD_ERROR_DECOMPOSITION.json':'MIXED_FIELD_ERROR',
              'P1_FIELD_PARTICLE_SENSITIVITY.json':'P1_FIELD_TRAJECTORY_IMPACT_LARGE'}
    rows={}
    for name,token in expected.items():
        data=read_json(ref/name);assert token in json.dumps(data),name
        rows[name]=dict(scientific_status=token,unchanged=True)
    scope=(root/'docs/original_scope/ACTIVE_DEVELOPMENT_SCOPE.md').read_text()
    assert 'STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT' in scope
    return dict(status='PASS_PRESERVED_FAILURES',reference_pipe=rows,P2P1='STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT',
                V3='PROPOSAL_ONLY',new_fixed_window_production='NOT_RUN')

def run_checks(output,root=ROOT,include_population=True):
    output=fresh_directory(output,root)
    checks={}
    tasks=[('preservation',preservation),('geometry',geometry.verify_geometry),('mesh',meshing.check_existing),('cfd_wss',cfd),('scientific_gates',scientific_gates)]
    if include_population:tasks.append(('population',population))
    for name,fn in tasks:
        print('START',name,flush=True)
        try:result=fn(root)
        except Exception as exc:
            result=dict(status='FAIL',error=repr(exc),traceback=traceback.format_exc())
        write_json(output/(name+'.json'),result);checks[name]=result['status'];print(name,result['status'],flush=True)
    write_json(output/'SUMMARY.json',dict(checks=checks,new_CFD_runs=0,new_formal_population_runs=0))
    return checks
