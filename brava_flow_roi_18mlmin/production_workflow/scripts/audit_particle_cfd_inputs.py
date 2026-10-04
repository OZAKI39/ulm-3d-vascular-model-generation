"""Read-only production input audit; does not import or integrate particles."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--identity', type=Path, required=True)
    args = ap.parse_args()
    identities = json.loads(args.identity.read_text())
    for name, row in identities.items():
        assert sha(args.source / name) == row['local_sha256'] == row['server_sha256'], name

    evidence = {}
    def cite(key, path, symbol, fragments):
        lines = (args.source / path).read_text().splitlines()
        hits = []
        for fragment in fragments:
            found = [i + 1 for i, line in enumerate(lines) if fragment in line]
            assert found, (path, fragment)
            hits.extend(dict(line=i, text=lines[i-1]) for i in found)
        evidence[key] = dict(module=path, function=symbol, sha256=sha(args.source/path), evidence=hits)
        return key

    cite('sample', 'field.py', 'FrozenFEMField.__init__/sample; FlowSample',
         ['g = np.einsum', 'self.curl_s_inv =', 'self.strain_s_inv =',
          'float(weights @ self.pressure_nodes_pa[nodes])', 'readonly(self.gradients_s_inv[cell])'])
    cite('cpu_free', 'population_v1/reference.py', 'ReferenceDynamics.provider/_trial_factory/advance',
         ['return sample.velocity_m_s, .5*sample.vorticity_s_inv',
          'gradient_at=lambda x:', 'return refine_interval'])
    cite('cpu_assembly', 'population_v1/components.py', 'assemble_components/solve_components',
         ['system=augment_planar_system', 'solved=(solve_resistance'])
    cite('planar', 'particle9a_motion.py', 'augment_planar_system',
         ['gap=gap_query(shape,wall)', 'gap.normal_inward,gradient_at(shape.center_m)', 'rim=rim_witness'])
    cite('shear', 'planar_wall_hydrodynamics.py', 'local_shear_vector/planar_wall_affine_block',
         ['return (np.eye(3)-np.outer(n,n))@g@n', 'force=SHEAR_FORCE_COEFFICIENT',
          'torque=SHEAR_TORQUE_COEFFICIENT', 'dr=L.T@dr4@L;db=L.T@db4'])
    cite('cpu_step', 'particle65_motion.py', 'v1_trial',
         ['free={i:np.concatenate(provider', 'handoff_constraints(old.shapes', 'end={i:s.moved(s.center_m+dt*v[i])',
          'wall_handoff_certificate(old.shapes[i]', 'wall_gaps[i]=wall_gap'])
    cite('refinement', 'physical_time_refinement.py', 'refine_interval',
         ['new,proof=trial(state,end-begin)', 'except TrialNeedsSubdivision', 'left,_=refine_interval'])
    cite('gpu_sample', 'population_v1/cuda/kernels.cu', 'fem_kernel',
         ['out[3]=', 'out[4+k]=gradient', 'out[13+k]=curl'])
    cite('gpu_step', 'population_v1/gpu_engine.py', 'GPUPopulationEngine._trial/_cpu_components',
         ['free=torch.cat([data[:,:3],.5*data[:,13:16]]', 'data[:,4:13].reshape(-1,3,3)',
          'velocity[ix],cpu_audit=self._cpu_components', 'end=p.x[:n]+dt*v', 'self.wall.path(p.x[:n]'])
    cite('gpu_resistance', 'population_v1/gpu_resistance.py', 'assemble_isolated/solve_isolated',
         ['normal=wall[:,7:10]', 'matrix=torch.diag_embed(diag);rhs=diag*free',
          'gn=(gradient@normal', 'shear=gn-normal', 'return velocity,hydro'])
    cite('gap', 'wall_gap.py', 'wall_gap; WallGapResult',
         ['gap_m: float', 'normal_inward: np.ndarray', 'nearest,_=wall.nearest_center_triangle',
          'capsule_triangle_many(capsule,wall.triangles[exact_ids])'])
    cite('nearfield', 'nearfield_regularization.py', 'NearFieldRegularizationV1.evaluate/lower_handoff_gap',
         ['h_lower_m=max(molecular,suspension)', 'eff=self.effective_gap(h', 'leading,n,old=evaluate_block'])
    cite('safety', 'population_v1/safety.py', 'verify_original_support',
         ['free={1:np.r_[sample.velocity_m_s,.5*sample.vorticity_s_inv]}',
          'lambda p:env.field.sample(p).velocity_gradient_s_inv'])

    def quantity(name, dtype, rhs, diagnostic, role, refs):
        return dict(name=name, data_type=dtype, directly_enters_dynamics=rhs,
                    diagnostic_or_transport_only=diagnostic, role=role, source_evidence=refs)
    quantities = [
        quantity('u', 'float64[3] m/s; CUDA float64[N,25][:,:3]', True, False,
                 'free translational velocity; Stokes diagonal forcing', ['sample','cpu_free','gpu_step','gpu_resistance']),
        quantity('G=du_i/dx_j', 'float64[3,3] 1/s; CUDA columns 4:13', True, False,
                 'piecewise constant raw P1 gradient; tangential shear s=(I-nn^T)G n', ['sample','planar','shear','gpu_step','gpu_resistance']),
        quantity('curl(u)', 'float64[3] 1/s; CUDA columns 13:16', True, False,
                 'derived from G; free angular velocity curl(u)/2 affects coupled resistance', ['sample','cpu_free','gpu_step']),
        quantity('WSS', 'not a FlowSample or GPU sampler field; separate CFD wall diagnostic', False, True,
                 'CFD_DERIVATIVE_DIAGNOSTIC; no WSS reader in production RHS', ['sample','gpu_sample','cpu_free','shear','gpu_step']),
        quantity('pressure', 'float Pa; CUDA column 3', False, True,
                 'sampled/transported for field API; no pressure or pressure-gradient force', ['sample','gpu_sample','cpu_free','gpu_step']),
        quantity('strain_rate', 'float64[3,3] 1/s; CUDA columns 16:25', False, True,
                 'sampled/transported; not separately read by spherical MB RHS', ['sample','gpu_sample','cpu_free','gpu_step']),
        quantity('wall gap', 'float m in WallGapResult; CUDA wall[:,0]', True, False,
                 'finite sphere-to-finite-triangle geometry; normal resistance, affine wall weight and handoff constraints', ['gap','nearfield','planar','gpu_resistance','cpu_step']),
        quantity('wall normal', 'float64[3] unit inward feature normal; CUDA wall[:,7:10]', True, False,
                 'exact nearest finite feature normal, local tangent frame, normal resistance and contact inequalities', ['gap','planar','shear','gpu_resistance']),
        quantity('wall curvature', 'no scalar/tensor curvature input', False, False,
                 'not explicitly evaluated; geometry changes closest feature/normal along a trajectory', ['gap','planar','gpu_resistance']),
        quantity('cell containment / cell ID', 'bool / int64', False, False,
                 'sampling and acceptance prerequisite; outside rejected, no extrapolation', ['sample','cpu_step','gpu_step']),
        quantity('wall feature/rim topology', 'triangle ID, feature category, rim flag, contact normals', True, False,
                 'selects local planar correction and exact constraints; rim fallback is declared production behavior', ['gap','planar','gpu_resistance']),
        quantity('mu / particle radius / pair geometry', 'float SI / float SI / gap+normal', True, False,
                 'material/source/geometry inputs, not derived from CFD solution; singleton benchmark has no pair partner', ['nearfield','planar','gpu_resistance']),
    ]
    probe = json.loads((ROOT/'manufactured_pipe/PARTICLE_PROBE_CONTRACT.json').read_text())
    cumulative = 0.
    for item in probe['distribution_contract']['bins']:
        mass = item['retained_probability']
        if cumulative + mass >= .5:
            median_um = .5*(item['low_um'] + (.5-cumulative)/mass*(item['high_um']-item['low_um']))
            median_bin = item['original_bin']
            break
        cumulative += mass
    result = dict(schema='PARTICLE_CFD_INPUT_CONTRACT_V1', status='AUDITED_BEFORE_NEARWALL_NUMERICAL_RESULTS',
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        method='Static call graph plus local/server byte identity; no trajectories/parity tests rerun',
        source_local_root=str(args.source.resolve()),
        source_server_root='/workspace/brava_flow_roi_18mlmin_20260928/population_v1/source/particle_3d/src/particle_3d',
        source_identity=identities, evidence=evidence, quantities=quantities,
        WSS_classification='CFD_DERIVATIVE_DIAGNOSTIC', WSS_trajectory_field_hard_gate=False,
        hard_gates=dict(particle_velocity_relative_L2=.05, full_G_relative_Frobenius_L2=.08,
            tangential_shear_vector_relative_L2=.08, vorticity_relative_L2=.08,
            bulk_velocity_relative_L2=.05, internal_flux_relative_error=.05, pressure_gradient_relative_error=.08),
        derivative_gate_definition='Each fixed particle band independently: full raw G, (I-nn^T)G n and curl(G); analytic cylinder n used for both fields to isolate field error. No per-zero-component percentages.',
        diagnostic_only_wall_distance_m=1e-7, hard_particle_wall_distances_m=probe['particle_relevant_wall_distances_m'],
        radius_contract=dict(distribution_sha256=probe['distribution_contract_sha256'],
            minimum_radius_m=probe['minimum_modeled_radius_m'], median_radius_m=median_um*1e-6,
            median_original_bin=median_bin, median_method='exact piecewise-uniform conditional inverse CDF at probability 0.5; no sampling'),
        wall_handoff=dict(h_lower_m=2e-9, formula='max(2e-9,0.001*a)', source_sha256=probe['nearfield_source_sha256']),
        integrator=dict(cpu='ReferenceDynamics -> p65.v1_trial -> refine_interval',
            gpu='GPUPopulationEngine._trial + event handling, binary subdivision and CPU component fallback',
            method='overdamped explicit accepted-time position update x_new=x+dt*v_RHS; same near-wall/contact rules',
            modifies_particle_physics=False, canonical_singleton_uses_production_engine=True),
        no_new_simulation=True)
    (ROOT/'PARTICLE_CFD_INPUT_CONTRACT.json').write_text(json.dumps(result, indent=2)+'\n')
    rows=[]
    for q in quantities:
        refs=[]
        for ref in q['source_evidence']:
            e=evidence[ref];n=e['evidence'][0]['line'];p=args.source/e['module']
            refs.append(f"[{e['module']}:{n}]({p.resolve()}:{n}) / `{e['function']}`")
        rows.append('| '+ ' | '.join([q['name'],q['data_type'], '是' if q['directly_enters_dynamics'] else '否', q['role'], '<br>'.join(refs)])+' |')
    report=f'''# Production particle CFD 输入审计

**审计已完成；在任何 NW 网格/CFD 结果出现前，固定 WSS=CFD_DERIVATIVE_DIAGNOSTIC。** WSS 仍报告，但不作为本轮 trajectory field hard gate。原 P0/P1 标定记录保持原结论，不追溯改成 PASS。这里没有运行新轨迹，也不证明时间步或物理精度。

已逐字节核对 {len(identities)} 个本地/服务器活动入口及核心文件，相同。范围为当前 BraVa `population_v1` CPU/GPU 生产路径与其显式阻力、近壁、采样依赖；不是其它 RBC 或历史 H0 模型。机器可读合同含每项源码 SHA、函数和精确行证据：[PARTICLE_CFD_INPUT_CONTRACT.json](PARTICLE_CFD_INPUT_CONTRACT.json)。

## 调用链

CPU：particle state → `ReferenceDynamics.provider` → `FrozenFEMField.sample` → `[u, curl(u)/2]` → `components.assemble_components` → 原 `assemble_v1` + `augment_planar_system` → `solve_resistance/solve_components` → `v1_trial` 的 `x+dt*v` → `refine_interval`。拒绝子步不消耗物理时间；接触事件及真实间隙约束保留。

GPU：pool state → `GPUFrozenFEMField.sample/fem_kernel` → FP64 u/G/curl → `assemble_isolated` → `solve_isolated` → `GPUPopulationEngine._trial`。困难组件回退到同一 CPU assembly；几何路径/终点检查后才接受位置和时间。`gpu_batch.prefetch` 复用同一装配与求解。`safety.verify_original_support` 复查也只读取 u/curl/G。这里的 integrator 是实际 engine 路径，不凭名叫 `integrator.py` 的旧模块推断。

| 输入 | 类型 | 直接进入动力学 | 用途/分类 | 源模块、函数与行 |
|---|---|---|---|---|
'''+ '\n'.join(rows)+f'''

## 为什么 WSS 不等于粒子读取的剪切

生产近壁仿射力读取球心的 `s=(I−nnᵀ)G n`，自由角速度读取 `curl(u)/2`；不读取壁面 `mu*(I−nnᵀ)*(G+Gᵀ)n` 或其模长。FlowSample 无 WSS 字段；CUDA 25 列依次是 u、p、G、curl、strain，而 RHS 只取 u/G/curl。压力与 strain 被采样不意味着进入力方程。法向和 gap 来自有限壁面几何，是额外的强依赖，不能由场精度通过推定几何动力学可靠。

## 预声明场门槛

每个冻结 particle band（0.377/1/2/5/10 µm）分别要求速度向量相对 L2≤5%，完整 G 的 Frobenius 相对 L2≤8%，剪切向量 s 相对 L2≤8%，curl 相对 L2≤8%。完整 G 是保守的总体检查，s 和 curl 显式覆盖实际使用的导数量；旧的单一轴向 du/dn 不足以替代它们。G_h 与 G_exact 用相同解析圆柱法向投影，避免把法向离散误差混入场门槛。零解析分量不逐项除零；整体向量/张量范数归一化，并保留绝对误差。bulk u≤5%，内部 Q≤5%，dp/dz≤8%。0.1 µm 和 WSS 单列诊断，探针不移动、不删除、不外推。

## 后续 canonical 模型身份

正式条件分布 SHA `{probe['distribution_contract_sha256']}`；a_min=0.375 µm，精确分段均匀逆 CDF 得 median radius={median_um:.12f} µm（原 bin {median_bin}），未抽样。两种半径的 h_lower 均为 2 nm。测试 gap 固定为 2 nm/0.5/2/5/10 µm。只允许替换背景场和 oracle 几何接口；粒子阻力、力律、有限尺寸、近壁、接触与积分自适应规则保持。当前未执行 oracle、CFD 轨迹或 dt 收敛；须先有合格场。
'''
    (ROOT/'PARTICLE_CFD_INPUT_AUDIT_ZH.md').write_text(report)
    print(json.dumps(dict(status=result['status'], source_files=len(identities), median_radius_um=median_um)))


if __name__ == '__main__':
    main()
