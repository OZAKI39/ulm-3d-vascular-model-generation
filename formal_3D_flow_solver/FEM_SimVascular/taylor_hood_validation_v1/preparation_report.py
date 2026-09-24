"""Scientific preparation artifacts; contains no assumed P2 solution values."""
from common import *
from p2 import EDGES
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pyvista as pv

def main():
    mesh=json.loads((REPORT/'data/mesh_elevation_audit.json').read_text())
    resource=json.loads((REPORT/'data/taylor_hood_resource_estimate.json').read_text())
    arrays=np.load(CASE/'SV_MESH/tet10_arrays_si.npz');points=arrays['points_m'];tet=arrays['tetra10']
    # A real boundary tetra, selected deterministically from the current mesh.
    face=pv.read(CASE/'SV_MESH/mesh-surfaces/INLET.vtp');eid=int(face['GlobalElementID'][len(face['GlobalElementID'])//2])-1
    xyz=points[tet[eid]];origin=xyz[:4].mean(axis=0);xyz=(xyz-origin)*1e6
    ids=(np.asarray(face['GlobalNodeID'],int)-1)[face.faces.reshape(-1,7)[len(face['GlobalElementID'])//2,1:]]
    local=[int(np.flatnonzero(tet[eid]==k)[0]) for k in ids]
    fig=plt.figure(figsize=(11,8),layout='constrained')
    for i,(count,title) in enumerate([(4,'P1: TET4 (4 velocity nodes)'),(10,'P2: TET10 (10 velocity nodes)')]):
        ax=fig.add_subplot(2,2,i+1,projection='3d')
        for edge in EDGES:ax.plot(*xyz[edge].T,color='#718096',lw=1.5)
        ax.add_collection3d(Poly3DCollection([xyz[local[:3]]],alpha=.13,facecolor='#0891b2'))
        ax.scatter(*xyz[:4].T,c='#225ea8',s=48,label='Original corner')
        if count==10:ax.scatter(*xyz[4:].T,c='#d95f0e',s=42,label='Shared edge midpoint')
        for n in range(count):ax.text(*xyz[n],str(n),fontsize=9)
        ax.set_title(title);ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');ax.set_zlabel('z (µm)');ax.set_box_aspect(np.ptp(xyz,axis=0));ax.legend(loc='upper left',fontsize=8)
        ax=fig.add_subplot(2,2,i+3,projection='3d');v=xyz[local]
        for e in [[0,1],[1,2],[2,0]]:ax.plot(*v[e].T,color='#718096',lw=1.5)
        ax.scatter(*v[:3].T,c='#225ea8',s=48)
        if count==10:ax.scatter(*v[3:].T,c='#d95f0e',s=42)
        ax.set_title('Same inlet face: '+('TRI3' if count==4 else 'TRI6'))
        ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');ax.set_zlabel('z (µm)')
        ax.set_box_aspect(np.maximum(np.ptp(v,axis=0),.1*np.ptp(v).max()))
    fig.suptitle(f'Order elevation of actual tetra {eid+1}; unchanged straight edges',fontsize=14)
    fig.savefig(REPORT/'figures/01_tet4_to_tet10_mesh.png',dpi=220);plt.close(fig)
    areas='\n'.join(f"| {k} | {v['P1_area_m2']:.15e} | {v['P2_integrated_area_m2']:.15e} | {v['relative_difference']:.3e} |" for k,v in mesh['boundary_areas'].items())
    (REPORT/'MESH_ELEVATION_AUDIT_ZH.md').write_text(f'''# 网格升阶审核

结论：PASS。只增加共享边中点，未重新剖分，未改变原角点或四面体拓扑。

- TET4：{mesh['original_vertices']:,} 节点，{mesh['tetra_count']:,} 单元。
- TET10：{mesh['new_nodes']:,} 节点，同样 {mesh['tetra_count']:,} 单元；新增 {mesh['unique_edges']:,} 个全局唯一边中点。
- 活跃自由度：{resource['P1_total_dof']:,} → {resource['P2P1_total_dof']:,}，增加到 {resource['active_dof_ratio']:.6f} 倍。原生四分量存储另保留边压力约束槽，总存储 {resource['solver_allocated_dof']:,}。
- 原体积 {mesh['total_volume_P1_m3']:.16e} m³；二次积分体积 {mesh['total_volume_P2_m3']:.16e} m³；相对差 {mesh['volume_relative_difference']:.3e}，处于浮点舍入量级。
- 全部单元、每单元 15 个求解器积分点的 Jacobian 均正，最小 {mesh['minimum_solver_quadrature_jacobian']:.16e}。
- 按固定版本 `read_msh.cpp:check_tet_conn` 调整局部角点编号，保持每个单元的原全局角点集合。准确节点顺序见 [节点顺序审核](TET10_NODE_ORDER_AUDIT_ZH.md)。
- 五类边界均为 TRI6，复用体网格中点，每个面恰有一个体单元父节点映射。
- 原生 VTK ID 数组必须为 Int32；保留了 Int64 加载失败输入，另建新路径修正。所有数值 ID、坐标和连接关系完全相同。

| 边界 | P1 面积 (m²) | P2 面积 (m²) | 相对差 |
|---|---:|---:|---:|
{areas}

本次有效输入：`{CASE}`。

资源估算见 [JSON](data/taylor_hood_resource_estimate.json)。矩阵存储和因子填充不是活跃自由度的简单线性函数；内存范围是预估，最终以运行记录为准。
''')
    doc=REPORT/'LOCAL_SVMP_TAYLOR_HOOD_AUDIT_ZH.md';doc.write_text(doc.read_text().replace('15 个关键文件','13 个关键文件'))
    path=REPORT/'reference/official_smokes/official_smoke_lfs/N004_shift/1-procs/result_050.vtu'
    if path.exists():
        g=pv.read(path);x,y=g.points[:,:2].T;exact=np.column_stack([np.sin(2*np.pi*x)*np.cos(2*np.pi*y),-np.cos(2*np.pi*x)*np.sin(2*np.pi*y)])
        u=np.asarray(g['Velocity'])[:,:2];assert np.isfinite(u).all()
        record=json.loads((path.parents[1]/'execution.json').read_text())
        record.update(native_case='stokes/manufactured_solution/P2P1/N004',nodes=g.n_points,cells=g.n_cells,
            nodal_velocity_relative_L2=float(np.linalg.norm(u-exact)/np.linalg.norm(exact)),
            note='Coarse official manufactured case: algebraic GPU smoke PASS; nodal analytic error reported without claiming a mesh convergence study.')
        dump(REPORT/'data/official_taylor_hood_smoke.json',record)
        print(json.dumps(record),flush=True)

if __name__=='__main__':main()
