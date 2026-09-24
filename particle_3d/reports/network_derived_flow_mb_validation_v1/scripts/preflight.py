"""Read-only frozen field/cohort audit; all writes stay in this new report."""
from pathlib import Path
import csv, hashlib, json, shutil, sys, time
import numpy as np
import pyvista as pv

R = Path(__file__).resolve().parents[1]
P = R.parents[1]
REPO = P.parent
F = Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular')
sys.path.insert(0, str(P/'src'))
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.wall_geometry import WallGeometry
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.inlet_flux import frozen_boundary_flux
from particle_3d.particle8_replay import canonical_hash
from particle_3d.wall_gap import wall_gap
from particle_3d.particle_shapes import Sphere

def sha(p):
    with Path(p).open('rb') as f: return hashlib.file_digest(f, 'sha256').hexdigest()
def dump(p, d):
    Path(p).write_text(json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
def compare(a,b):
    a,b=np.asarray(a),np.asarray(b)
    return dict(shape_old=list(a.shape),shape_new=list(b.shape),dtype_old=str(a.dtype),dtype_new=str(b.dtype),
        bitwise_equal=a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes(),
        numeric_equal=np.array_equal(a,b),sha_old=hashlib.sha256(a.tobytes()).hexdigest(),sha_new=hashlib.sha256(b.tobytes()).hexdigest())
def write_csv(path, rows):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    start=time.time(); oldcase=F/'flow_cases/mean-2p0-mmps';newcase=F/'flow_cases/mean-2p0-mmps-A-H0-pressure-v1'
    approval=newcase/'reports/physics_validation_H0.json'; approved=json.loads(approval.read_text())
    assert approved['status']=='PASS' and approved['export']['status']=='PASS'
    names=[n for n in approved['export']['files'] if n.endswith('.vtu')]; assert len(names)==1
    paths={'OLD':oldcase/'frozen_flow/steady_flow_mean_2p0_mmps.vtu','NEW':newcase/'frozen_flow'/names[0]}
    assert sha(paths['NEW'])==approved['export']['files'][names[0]]['sha256']
    grids={k:pv.read(p) for k,p in paths.items()};a,b=grids.values()
    contract=dict(authoritative_new_source=str(approval),authoritative_new_source_sha256=sha(approval),
        paths={k:str(v) for k,v in paths.items()},sha256={k:sha(v) for k,v in paths.items()},
        point_count={'OLD':a.n_points,'NEW':b.n_points},cell_count={'OLD':a.n_cells,'NEW':b.n_cells},
        coordinates=compare(a.points,b.points),tetra_connectivity=compare(a.cells,b.cells),cell_types=compare(a.celltypes,b.celltypes),boundaries={})
    for role in ['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        sp=[case/'SV_MESH/mesh-surfaces'/f'{role}.vtp' for case in [oldcase,newcase]];s,t=map(pv.read,sp)
        checks={'coordinates':compare(s.points,t.points),'faces':compare(s.faces,t.faces)}
        for assoc in ['point_data','cell_data']:
            assert set(getattr(s,assoc))==set(getattr(t,assoc))
            for key in getattr(s,assoc):checks[assoc+'/'+key]=compare(getattr(s,assoc)[key],getattr(t,assoc)[key])
        contract['boundaries'][role]=dict(paths=list(map(str,sp)),file_sha256=list(map(sha,sp)),checks=checks)
    checks=[contract[n] for n in ['coordinates','tetra_connectivity','cell_types']]+[c for row in contract['boundaries'].values() for c in row['checks'].values()]
    contract['geometry_identical']=all(c['numeric_equal'] for c in checks)
    contract['geometry_bitwise_identical']=all(c['bitwise_equal'] for c in checks)
    contract['face_ids']={'WALL':1,'OUTLET_03':2,'OUTLET_01':3,'INLET':4,'OUTLET_02':5}
    contract['status']='PASS' if contract['geometry_identical'] else 'OLD_NEW_FLOW_GEOMETRY_MISMATCH'
    dump(R/'data/old_new_flow_contract.json',contract)
    if not contract['geometry_identical']: raise RuntimeError(contract['status'])
    provenance,mesh,historical_flow,boundaries=read_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    assert provenance['flow_sha256']==contract['sha256']['OLD']
    for role, s in boundaries.items():
        t=pv.read(oldcase/'SV_MESH/mesh-surfaces'/f'{role}.vtp')
        assert np.array_equal(s.points,t.points) and np.array_equal(s.faces,t.faces)
        for name in s.point_data:assert np.array_equal(s.point_data[name],t.point_data[name])
    fields={k:FrozenFEMField.from_grids(mesh,v) for k,v in grids.items()}
    flux={}
    for name,g in grids.items():
        audit,_=frozen_boundary_flux(mesh,g,boundaries)
        qi=audit['INLET']['signed_Q_m3_s']
        flux[name]=dict(integration='Exact P1 nodal triangle integration of authoritative VTU',boundaries=audit,Qin_m3_s=qi,
            split={f'O{i}':audit[f'OUTLET_{i:02d}']['signed_Q_m3_s']/qi for i in range(1,4)},
            global_mass_balance_pass=abs(audit['mass_balance']['relative_signed_residual'])<1e-6)
        assert flux[name]['global_mass_balance_pass']
    dump(R/'data/reintegrated_flow_flux.json',flux)
    wall=WallGeometry.from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    history=P/'outputs/particle9a1_2mmps'; configpath=history/'provenance/smoke30_config.json';ledgerpath=history/'admission/birth_ledger.json'
    config=json.loads(configpath.read_text());ledger=json.loads(ledgerpath.read_text());events=ledger['events'][:30]
    assert len(events)==30 and set(map(str,[e['particle_id'] for e in events]))==set(config['event_sha256'])
    source_mismatch=[n for n,h in config['p9a1_identity']['source_sha256'].items() if sha(P/'src/particle_3d'/n)!=h]
    assert not source_mismatch
    assert config['p9a1_identity']['flow']['flow_sha256']==contract['sha256']['OLD']
    rows=[];admissibility=[];initial={};copied=[];out=R/'outputs/OLD/trajectories';out.mkdir(parents=True)
    for e in events:
        pid=e['particle_id'];assert e['admission_strategy']=='B'
        assert canonical_hash(e)==config['event_sha256'][str(pid)]
        stem=history/'P9A1/trajectories'/f'mb_{pid:06d}';meta=json.loads(stem.with_suffix('.json').read_text())
        assert meta['birth_metadata']==e and meta['birth_metadata_sha256']==canonical_hash(e)
        assert meta['p9a1_identity']==config['p9a1_identity']
        assert meta['samples_sha256']==sha(stem.with_suffix('.npz'))
        samples=np.load(stem.with_suffix('.npz'))['samples'];assert np.array_equal(samples[0,1:4],e['birth_center_m'])
        assert np.array_equal(samples[0,10:14],e['q'])
        initial[str(pid)]=dict(velocity_m_s=samples[0,4:7].tolist(),omega_s_inv=samples[0,7:10].tolist(),q=samples[0,10:14].tolist(),position_m=samples[0,1:4].tolist())
        pair={}
        for label,field in fields.items():
            particle,status,detail=FiniteSizeAdmission(None,None,wall=wall,field=field).check(e,np.array(e['birth_center_m']),{})
            sample=field.sample(e['birth_center_m']);gap=wall_gap(Sphere(e['birth_center_m'],e['radius_m']),wall,inside_lumen=True)
            pair[label]=dict(status=status,accepted=particle is not None,inside_fluid=bool(sample.inside_lumen),initial_tetra=sample.tetra_id,
                wall_gap_m=gap.gap_m,geometric_detail=detail,cap='INLET',anchor_triangle=e['anchor_triangle'],open_boundary_rule='WALL only is solid; original P9-A.1 rim topology',
                field_initial_diagnostic_velocity_m_s=sample.velocity_m_s.tolist(),field_initial_diagnostic_omega_s_inv=(.5*sample.vorticity_s_inv).tolist())
        equal=all(pair['OLD'][k]==pair['NEW'][k] for k in ['status','accepted','inside_fluid','initial_tetra','wall_gap_m','geometric_detail','cap','anchor_triangle','open_boundary_rule'])
        assert equal and pair['OLD']['accepted'],'FLOW_FIELD_GEOMETRY_CONTRACT_VIOLATION'
        admissibility.append(dict(bubble_id=pid,identical=equal,**pair))
        row=dict(bubble_id=pid,diameter_um=e['diameter_um'],radius_m=e['radius_m'],radius_um=e['radius_m']*1e6,
            initial_x_m=e['birth_center_m'][0],initial_y_m=e['birth_center_m'][1],initial_z_m=e['birth_center_m'][2],initial_tetra=pair['OLD']['initial_tetra'],
            initial_velocity_m_s=json.dumps(initial[str(pid)]['velocity_m_s']),initial_omega_s_inv=json.dumps(initial[str(pid)]['omega_s_inv']),
            initial_q=json.dumps(e['q']),position_seed=json.dumps(e['position_seed']),anchor_triangle=e['anchor_triangle'],
            inlet_sample_quantile='NOT_STORED_HISTORICALLY; exact anchor and seed recovered',event_sha256=canonical_hash(e),all_deterministic_metadata_json=json.dumps(e,sort_keys=True))
        rows.append(row)
        for suffix in ['.json','.npz','.audit.jsonl.gz']:
            src=stem.with_suffix(suffix);assert src.is_file();dest=out/src.name;shutil.copy2(src,dest)
            copied.append(dict(source=str(src),copy=str(dest),sha256=sha(src)))
    write_csv(R/'data/paired_cohort_manifest.csv',rows)
    dump(R/'data/paired_events.json',events);dump(R/'data/paired_initial_states.json',initial)
    dump(R/'data/cohort_geometry_admissibility.json',dict(status='PASS',N=30,wall_roundoff_m=wall.roundoff_m,rows=admissibility))
    dump(R/'data/cohort_provenance.json',dict(source='HISTORICAL_P9A1_SMOKE30_EXACT_REUSE',N=30,master_seed=ledger['master_seed'],
        source_ledger=str(ledgerpath),source_ledger_sha256=sha(ledgerpath),source_config=str(configpath),source_config_sha256=sha(configpath),
        historical_source_mismatches=source_mismatch,historical_science_files_verified=len(config['p9a1_identity']['source_sha256']),
        OLD_reused=True,OLD_files=copied,manifest_sha256=sha(R/'data/paired_cohort_manifest.csv'),events_sha256=sha(R/'data/paired_events.json'),
        initial_states_sha256=sha(R/'data/paired_initial_states.json'),admission_method='Historical Method B external inlet; no new generation, no Method C or interior injection',
        initial_velocity_policy='Freeze saved OLD row-zero velocity/omega/q in both sets; overdamped subsequent steps use the respective background field',
        missing_historical_rng_quantile='Quantile not saved; exact anchor, event, stream seed and full metadata preserved, no resampling',elapsed_seconds=time.time()-start))
    print(json.dumps(dict(geometry_identical=True,flux=flux,N=30,old_reused=True,elapsed_seconds=time.time()-start)),flush=True)
if __name__=='__main__':main()
