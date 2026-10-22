"""Raw P1 u/G/shear/curl at the immutable manufactured probes; no CFD launch."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent/'vendor'))
from pipe_measurements import PipeMeasurements
from pipe_probe_geometry import frozen_probes, strict_cells
from poiseuille import reference


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def curl(g):
    return np.stack([g[:,2,1]-g[:,1,2],g[:,0,2]-g[:,2,0],g[:,1,0]-g[:,0,1]],axis=1)


def shear(g,n):
    gn=np.einsum('nij,nj->ni',g,n)
    return gn-n*np.sum(n*gn,axis=1)[:,None]


def errors(actual,exact):
    delta=(actual-exact).reshape(len(actual),-1)
    norm=float(np.linalg.norm(delta));den=float(np.linalg.norm(exact))
    return dict(relative_L2=norm/den if den else None,absolute_L2=norm,
                absolute_RMS=norm/np.sqrt(len(actual)),
                max_absolute_error=float(np.linalg.norm(delta,axis=1).max()),
                reference_L2=den)


def analyze(name):
    case=ROOT/'manufactured_pipe'/name/'cfd'
    inp=ROOT/'straight_pipe/manufactured'/name/'INPUT.json'
    inputs=json.loads(inp.read_text());proto=json.loads((ROOT/'manufactured_pipe/PROTOCOL.json').read_text())
    execution=json.loads((case/'reports/execution.json').read_text());assert execution['status']=='PASS'
    previous=json.loads((case/'analysis/CALIBRATION.json').read_text())
    field_path=case/execution['final_vtu'];mesh_path=case/'SV_MESH/mesh-complete.mesh.vtu'
    for p in [field_path,mesh_path]:assert sha(p)==previous['source_sha256'][str(p.relative_to(ROOT))]
    m=PipeMeasurements(case/'SV_MESH/mesh_arrays.npz',inputs['analytic']['Q_m3_s'],proto['mean_velocity_m_s'])
    u,p=m.read(field_path);points,tets=m.points,m.tetra
    dx=points[tets[:,1:]]-points[tets[:,:1]]
    G=np.linalg.solve(dx,u[tets[:,1:]]-u[tets[:,:1]]).transpose(0,2,1)
    levels,x,groups,planes,angles=frozen_probes(proto,inputs['radius_m'],inputs['length_m'])
    frozen=np.load(case/'analysis/FIXED_PROBES.npz')
    assert np.array_equal(x,frozen['points_m'])
    cells,locators=strict_cells(pv.read(mesh_path),points,tets,x)
    assert np.array_equal(cells,frozen['cell']), 'Frozen containment decision changed'
    assert np.array_equal(G,frozen['gradient_all_tetra']), 'Raw P1 gradient changed'
    valid=cells>=0;uh=np.full((len(x),3),np.nan);gh=np.full((len(x),3,3),np.nan)
    c=cells[valid]
    b=np.linalg.solve(dx[c].transpose(0,2,1),(x[valid]-points[tets[c,0]])[...,None])[...,0]
    weights=np.c_[1-b.sum(axis=1),b]
    uh[valid]=np.einsum('ni,nij->nj',weights,u[tets[c]]);gh[valid]=G[c]
    assert np.array_equal(uh,frozen['velocity'])
    n=-np.c_[np.cos(angles),np.sin(angles),np.zeros(len(x))]
    exact=reference(x,inputs['radius_m'],inputs['length_m'],proto['mean_velocity_m_s'],proto['mu_Pa_s'])
    actual=dict(u=uh,G=gh,s=shear(gh,n),curl=curl(gh))
    expected=dict(u=exact['u_m_s'],G=exact['G_s_inv'],s=shear(exact['G_s_inv'],n),curl=curl(exact['G_s_inv']))
    rows=[]
    for group,d in enumerate(levels):
        mask=groups==group;role=proto['probe_roles'][group];complete=bool(valid[mask].all())
        metrics={k:errors(actual[k][mask],expected[k][mask]) if complete else None for k in actual}
        old=previous['samples'][group]
        if complete:assert abs(metrics['u']['relative_L2']-old['velocity_relative_L2'])<1e-14
        rows.append(dict(wall_distance_m=d,role=role,samples=int(mask.sum()),inside=int(valid[mask].sum()),
            complete_coverage=complete,metrics=metrics,old_du_dn_relative_L2=old['gradient_relative_L2'],
            WSS_diagnostic_relative_L2=old['raw_wall_WSS_relative_L2'],
            input_gates={k:complete and metrics[k]['relative_L2']<=limit for k,limit in [('u',.05),('G',.08),('s',.08),('curl',.08)]} if role=='PARTICLE_RELEVANT' else None))
    return dict(case=name,status='LEGACY_POSTPROCESS_ONLY',historical_status_unchanged=previous['status'],
        current_particle_input_gates_pass=all(all(q['input_gates'].values()) for q in rows if q['input_gates'] is not None),
        groups=rows,raw_WSS_relative_L2=previous['raw_WSS_relative_L2'],
        internal_sections=previous['interior_sections'],dp_dz_relative_error=previous['pressure_gradient_relative_error'],
        input_sha256={str(q.relative_to(ROOT)):sha(q) for q in [field_path,mesh_path,inp,case/'analysis/FIXED_PROBES.npz',case/'analysis/CALIBRATION.json']},
        probes_sha256=hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest(),
        raw_fields_and_saved_gradients_match=True,raw_velocity_and_cell_ids_match=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--legacy',action='store_true');args=ap.parse_args()
    assert args.legacy
    result=dict(schema='LEGACY_PARTICLE_INPUT_BASELINE_V1',CFD_runs=0,
        derivative_definition='G_ij=du_i/dx_j; inward analytic n=(-cos(theta),-sin(theta),0); s=(I-nn^T)G n; curl from same G',
        error_definition='absolute_L2=sqrt(sum over probes of squared vector/Frobenius error); RMS=L2/sqrt(N); max_absolute_error=max per-probe vector/Frobenius error',
        units=dict(u='m/s',G='1/s',s='1/s',curl='1/s'),
        smoothing=False,extrapolation=False,probe_changes=False,WSS_role='CFD_DERIVATIVE_DIAGNOSTIC',
        input_contract_sha256=sha(ROOT/'PARTICLE_CFD_INPUT_CONTRACT.json'),
        sources_sha256={n:sha(ROOT/'scripts'/n) for n in ['particle_field_inputs.py','pipe_probe_geometry.py','poiseuille.py','pipe_measurements.py']},
        cases=[analyze(n) for n in ['R1_P0','R1_P1']])
    (ROOT/'LEGACY_P0_P1_PARTICLE_INPUT_BASELINE.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    block=['\n<!-- LEGACY_INPUT_BASELINE_START -->','## 旧 P0/P1 的真实粒子输入基线（仅新增后处理）','',
        '直接读取已收敛原始场，重新严格定位并证明 probe 坐标、cell ID、u 和全部原始 P1 G 与保存数据一致。没有重新运行 CFD，旧 CALIBRATION.json 及其 FAIL 判定不变。下表每项均为相对 L2；五个 band 每组96点，G为完整张量，s/curl为完整向量。绝对 L2、RMS、最大逐点绝对误差及全部冻结探针见[基线 JSON](LEGACY_P0_P1_PARTICLE_INPUT_BASELINE.json)。','',
        '| Case | 壁距 µm | 覆盖 | u | G | s | curl | 旧 du/dn | WSS（诊断） |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for case in result['cases']:
        for q in case['groups']:
            if q['role']!='PARTICLE_RELEVANT':continue
            values=[q['metrics'][k]['relative_L2'] for k in ['u','G','s','curl']]+[q['old_du_dn_relative_L2'],q['WSS_diagnostic_relative_L2']]
            block.append(f"| {case['case']} | {q['wall_distance_m']*1e6:g} | {q['inside']}/{q['samples']} | "+' | '.join(f'{v*100:.4f}%' for v in values)+' |')
    block.extend(['','WSS 的分类和 u/G/s/curl 门槛均沿用已冻结合同；没有根据这些新数值改变门槛。','<!-- LEGACY_INPUT_BASELINE_END -->',''])
    path=ROOT/'PARTICLE_CFD_INPUT_AUDIT_ZH.md';text=path.read_text()
    if '<!-- LEGACY_INPUT_BASELINE_START -->' in text:text=text[:text.index('<!-- LEGACY_INPUT_BASELINE_START -->')].rstrip()+'\n'
    path.write_text(text+'\n'.join(block))
    print('\n'.join(block))


if __name__=='__main__':main()
