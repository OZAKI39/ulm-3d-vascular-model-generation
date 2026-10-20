"""Read accepted results, audit delivery, and write the scientific application report."""
from common import *
from pathlib import Path
import csv,time,datetime

def main():
    case=O/'cases/biphasic';summary=json.loads((case/'SUMMARY.json').read_text());freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
    for p,want in freeze['implementation_sha256'].items():assert sha(R/p)==want,p
    for p,want in freeze['input_sha256'].items():assert sha(O/p)==want,p
    ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);N=len(ref['X']);bt=ref['bond_types'];b=ref['bonds']
    field=load_npz(O/'BIPHASIC_FIELD.npz');source=load_npz(Y_SOURCE/'provenance/reference.npz')
    unchanged={k:bool(np.array_equal(ref[k],source[k])) for k in ['X','bonds','bond_types','original_bond_ids','r0','taper','m','V','initial_cut']};assert all(unchanged.values())
    graph=module('finish_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
    final=case/'keys'/('final.npz' if summary['status']=='COMPLETED' else 'last_accepted.npz');state=load_npz(final)
    lab,sizes,attached=graph.topology(ref,state['d']);free=sizes[~attached]
    assert np.array_equal(state['atom_ids'],np.arange(1,N+1)) and np.array_equal(state['x'][n:],ref['X'][n:])
    assert sizes.sum()==n and not state['d'][bt==3].any()
    with (case/'history.csv').open() as f:rows=list(csv.DictReader(f))
    def a(k):return np.array([float(z[k]) for z in rows])
    material=json.loads((O/'BIPHASIC_MATERIAL_SPEC.json').read_text());visual=json.loads((O/'VISUALIZATION_VALIDATION.json').read_text())
    assert visual['status']=='PASS'
    with (F/'FINAL_FRAGMENT_DISTRIBUTION.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['label','nodes','attached','volume_m3','mass_kg','equivalent_volume_diameter_um'])
        for i in range(len(sizes)):w.writerow([i,sizes[i],bool(attached[i]),sizes[i]*V,sizes[i]*RHO*V,(6*sizes[i]*V/np.pi)**(1/3)*1e6])
    bins={f'{lo}-{hi}':dict(components=int(((free>=lo)&(free<=hi)).sum()),nodes=int(free[(free>=lo)&(free<=hi)].sum())) for lo,hi in [(1,2),(3,19),(20,200),(201,1927)]}
    actual=dict(status=summary['status'],accepted_time_ms=summary['time_s']*1e3,material_nodes=n,anchor_nodes=N-n,
        detached_nodes=int(free.sum()),detached_fraction=float(free.sum()/n),free_components=len(free),largest_free_nodes=int(free.max(initial=0)),
        final_20_200=int(((free>=20)&(free<=200)).sum()),maximum_20_200=summary['max_medium_components'],fragment_bins=bins,
        failed_bulk=int(((state['d']>=1)&(bt==1)).sum()),failed_interface=int(((state['d']>=1)&(bt==2)).sum()),
        accepted_steps=summary['accepted_steps'],predictive_rejections=summary['predictive_rejections'],
        max_relative_energy_residual=summary['max_relative_energy_residual'],minimum_wall_clearance_m=float(a('min_wall_clearance_m').min()),
        maximum_wall_contact_nodes=int(a('wall_contact_nodes').max()),maximum_L2_patches=int(a('L2_active_patches').max()),
        final_W_drive_J=float(a('W_drive_J')[-1]),final_W_drag_J=float(a('W_drag_J')[-1]),final_W_wall_J=float(a('W_wall_J')[-1]),
        peak_actual_surface_traction_Pa=float(a('shear_peak_Pa').max()),final_wall_potential_J=float(state['wall_potential_J']),
        wall_s=summary['wall_s'],cpu_s=summary['cpu_s'],peak_RSS_bytes=summary['peak_RSS_bytes'],
        initial_geometry_and_material_graph_unchanged=unchanged,material_phase_fraction=material['strong_fraction'],
        original_node_membership_complete=True,native_library_sha256=summary['library_sha256'],stop_reasons=summary['stop_reasons'])
    dump(O/'APPLICATION_RESULT.json',actual)
    numerical_pass=summary['status']=='COMPLETED' and actual['max_relative_energy_residual']<=.01
    wall_clearance_ok=actual['minimum_wall_clearance_m']>=-2e-6
    dump(O/'DELIVERY_VALIDATION.json',dict(status='PASS' if numerical_pass and wall_clearance_ok else 'LIMITED_RESULT',
        numerical_protocol_completed=numerical_pass,wall_center_penetration_within_2um=wall_clearance_ok,
        frozen_sources_unchanged=True,original_Y_geometry_exact=True,mass_and_atom_ids_closed=True,visualization=visual['status'],
        verified_full_CFD_coupling=False,calibrated_material_microstructure=False))
    report=f'''# Y_CLOT_CURRENT_METHOD_APPLICATION_V1

运行 `{RUN}`。按用户要求，将当前 V3 表面耦合、强弱材料场和变分 LPS 软化用于原 Y 血管内不规则血块；仅执行一次正式工况。

## 实际状态

- 原生状态：**{summary['status']}**；接受物理时间 **{actual['accepted_time_ms']:.6f} ms**，计划 3.6 ms。
- 末态游离 **{actual['detached_nodes']}/{n} 点（{actual['detached_fraction']:.3%}）**，形成 **{actual['free_components']}** 个真实材料连通分量；最大 **{actual['largest_free_nodes']}** 点。
- 20–200 点游离块末态 **{actual['final_20_200']}** 个，历史最多 **{actual['maximum_20_200']}** 个。该区间只是统计，不是预切块或材料区设计目标。
- 本体/界面完全失效键：**{actual['failed_bulk']} / {actual['failed_interface']}**；固定支架没有损伤。

## 保留与适配

原 1927 血栓点、1561 锚点、157637 条原始材料键及全部参考 m/r0/V/taper 和 ID 逐项相同。原几何没有初始切断键；没有移入规则试样缺口，也没有按当前位置重新连键。此几何是历史工程构造的三鼓包贴壁血块，不是患者影像分割。

G=940 Pa、K=100000/3 Pa、ρ=1000 kg/m³、dx=50 µm、horizon=150.5 µm。弱/混合/强的 (s0,sc) 分别为 (0.10,0.12)/(0.175,0.210)/(0.25,0.30)，界面为 (0.140625,0.16875)。原生软化库复用，不修改本构或生产库。

材料场只做一次平移 (-6.0,0.05,1.0) mm 后按原血块位置采样，没有缩放、随机种子搜索或结果驱动调参。保留 8 个原强区标签、882 强点，比例 **45.7706%**；原规则试样为 33.3292%。裁剪导致比例和强区体积改变，不能把本轮与规则试样的差异归因于几何单一因素。强区物理尺度仍未标定。

L1 将规则格点的四面体支撑构造适配到缺点的不规则域，只有四个原节点和六条原始材料边均存在才保留辅助四面体。原血块图已排除跨越曲面空隙的几何近邻，不能补键来凑齐单元。初始化曾因直接查找不存在的原始边而失败，修正为上述完整支撑检查；这是正式启动前的适配错误，失败目录保留。L2 沿用原对偶控制体积、存活非局部支撑、开口与遮挡规则。

新增管腔检查要求当前面中心及外向 3 µm 探针位于原 Y 流体四面体域内。它不等于求解实际绕流或裂缝流动。初始通过 1590 个 L1 面，求积面积约 2.22567 mm²；完整连续湿润表面未认证。

## 载荷、血管壁与数值方法

沿用 +X 当前切平面投影、β=μ/h=40 Pa·s/m、名义 1000 Pa，0–2 ms 平滑升、2–2.4 ms 保持、2.4–2.8 ms 降、2.8–3.6 ms 被动相对阻力。载荷是机械代理；没有叠加旧 Gaussian、静态 CFD 节点力、COM 拖曳或数值阻尼，没有模拟 PZT 或微泡声流。

原 Y 壁几何进入接触：有效半径20 µm、刚度50 N/m、阻尼比0.5。弹性法向使用真实最近距离的梯度，使其与 1/2 k overlap² 势能一致；壁力与壁功单列。当前未做长时出口交接；若到达出口，本实现明确停止而非将出口当壁。最近轮 L1/L2 仍只是部分近似表面加载。

原生 Velocity-Verlet，基础0.25 µs、最小0.03125 µs，单步损伤增量≤0.1，能量余额≤1%。速度相关载荷在原生后半步内迭代一致求值；只有接受态提交表面历史。原始存活材料图判定附着和游离，接触不恢复材料连接。

| 验证/成本 | 实测 |
|---|---:|
| 接受步 / 预测细分拒绝 | {actual['accepted_steps']} / {actual['predictive_rejections']} |
| 最大归一化能量余额 | {actual['max_relative_energy_residual']:.6%} |
| 最小壁净距 | {actual['minimum_wall_clearance_m']*1e6:.6f} µm |
| 最大壁接触节点数 | {actual['maximum_wall_contact_nodes']} |
| 历史最大有效 L2 面数 | {actual['maximum_L2_patches']} |
| 最大实际表面牵引 | {actual['peak_actual_surface_traction_Pa']:.6f} Pa |
| 驱动功 / 流体阻力功 / 壁功 | {actual['final_W_drive_J']:.9e} / {actual['final_W_drag_J']:.9e} / {actual['final_W_wall_J']:.9e} J |
| CPU / 墙钟 | {actual['cpu_s']:.3f} / {actual['wall_s']:.3f} s |
| 峰值 RSS | {actual['peak_RSS_bytes']/1024**3:.4f} GiB |

能量账本为 Wdrive+Wdrag+Wwall=ΔU+ΔK+D+ΔUc−Jc+R。壁作用作为外功进入；不得再把壁势能重复加入右端。残差R不并入损伤耗散。数值准入不替代材料强度、Gc、界面黏附、局部牵引及完整水动力几何的实验验证。

## 可视化与复现

[打开图像与视频](../../../results/{TASK}/{RUN}/OPEN_RESULTS.html)。视频1920×1080、24fps，红色血栓、Z轴竖直、真实位移×1。所有画面来自接受态，未制造中间动力学。全景视频仅显示存活短键以降低遮挡，完整存活图仍用于碎块定义。受载面和箭头来自保存载荷；2µN箭头对应200µm图上世界长度。载荷图与CSV使用 `surface_force=surface_traction×geometric_area`，开放权重已在牵引中，不能重复乘入。

源码：`analysis/{TASK}/`；冻结输入和原生状态：`verification/{TASK}/{RUN}/`。正式运行由 supervisor 管理，不自动重启。初始微步和重启检查是接口检查，不是额外完整工况。使用 `SOURCE_FREEZE.json` 检查源码、输入和库身份；新复现必须新目录，不覆盖本次接受证据。

本次没有重新求解CFD，没有双向FSI，没有声学能力标定或AI控制。停止原因：{summary['stop_reasons']}。
'''
    (O/'REVIEW_Y_CLOT_CURRENT_METHOD_V1.md').write_text(report)
    (F/'DEVELOPMENT_REPORT.md').write_text(report.replace(f'../../../results/{TASK}/{RUN}/OPEN_RESULTS.html','OPEN_RESULTS.html'))
    dump(O/'FINAL_STATUS.json',dict(status=summary['status'],task=TASK,run=RUN,formal_runs=1,results=actual,
        report='REVIEW_Y_CLOT_CURRENT_METHOD_V1.md',figures=str(F/'OPEN_RESULTS.html'),visualization_pass=True,
        scientific_label='EXPLORATORY_TRANSFERRED_Y_CLOT_PRESCRIBED_SHEAR',CFD_recomputed=False,GitHub_pushed=False))
    print('Y_APPLICATION_REPORT_COMPLETE',actual,flush=True)

if __name__=='__main__':main()
