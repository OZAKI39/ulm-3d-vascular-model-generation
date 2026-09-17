"""
把 SWC 树转换成圆柱面片 mesh。

SWC 只保存点和半径；三维显示需要把每条父子边变成一个圆柱。这里直接
构造 NumPy 顶点数组和 face 数组，最后交给 PyVista 生成 PolyData。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .swc_loader import SWCTree

RadiusStyle = Literal["cylindrical", "tapered"]
CapStyle = Literal["none", "terminal", "all"]


@dataclass(frozen=True)
class MeshData:
    """PyVista PolyData 所需的 NumPy mesh 缓冲区。"""

    # 所有顶点坐标，形状为 N x 3。
    points: np.ndarray
    # PyVista 面片数组；每个面以前缀数字表示顶点个数，例如 [4, a, b, c, d]。
    faces: np.ndarray
    # 每个面片对应的血管半径，用于按半径着色。
    cell_radii_um: np.ndarray
    # 每个面片来自哪条 SWC 边，便于调试和后处理。
    cell_edge_ids: np.ndarray
    # 原始 SWC 边总数，不含被跳过的退化边。
    edge_count: int
    # 长度接近 0 的边无法形成圆柱，会计入 skipped_edges。
    skipped_edges: int


def _basis_from_direction(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """给定圆柱轴向，构造两个垂直于轴向的单位基向量。"""
    # 先把方向向量归一化，后续叉乘只关心方向不关心长度。
    unit = direction / np.linalg.norm(direction)
    # 从 x/y/z 三个坐标轴里选一个最不平行于圆柱轴的参考轴。
    candidates = np.eye(3)
    # 点积绝对值越小，说明候选轴与圆柱轴越接近垂直，叉乘越稳定。
    reference = candidates[np.argmin(np.abs(candidates @ unit))]
    # 第一个截面基向量 u 与圆柱轴、参考轴都垂直。
    u = np.cross(unit, reference)
    # 归一化 u，确保半径乘以 u 后长度仍等于半径。
    u /= np.linalg.norm(u)
    # 第二个截面基向量 v 与圆柱轴和 u 垂直，组成圆截面平面。
    v = np.cross(unit, u)
    # 返回两个截面方向，后面用 cos/sin 组合出圆环。
    return u, v


def _append_cap(
    faces: list[int],
    cell_radii: list[float],
    cell_edge_ids: list[int],
    center_index: int,
    ring: list[int],
    radius_um: float,
    edge_id: int,
    reverse: bool,
) -> None:
    """给圆柱端面追加三角形 cap。"""
    # ring 是端面圆环上的顶点编号列表，n 就是圆周离散边数。
    n = len(ring)
    # 每两个相邻圆环点与中心点组成一个三角形，拼满整个端面。
    for i in range(n):
        # 当前圆环点。
        a = ring[i]
        # 下一个圆环点，最后一个点通过取模回到第一个点。
        b = ring[(i + 1) % n]
        # reverse 用于控制三角形顶点顺序，从而保持法向朝外。
        if reverse:
            faces.extend([3, center_index, b, a])
        else:
            faces.extend([3, center_index, a, b])
        # 端盖三角形也记录对应半径，方便颜色映射。
        cell_radii.append(radius_um)
        # 端盖属于当前这条血管边。
        cell_edge_ids.append(edge_id)


def build_cylinder_mesh(
    tree: SWCTree,
    *,
    radius_style: RadiusStyle = "cylindrical",
    sides: int = 18,
    radius_scale: float = 1.0,
    cap_style: CapStyle = "terminal",
) -> MeshData:
    """
    为整棵 SWC 树构建批量圆柱 mesh。

    `cylindrical`：父子节点之间用同一半径，匹配本生成器“一段 vessel 一个
    半径”的模型。
    `tapered`：父端用父节点半径，子端用子节点半径，更接近 Vaa3D 这类
    SWC 浏览器的显示方式。
    """
    # 先校验枚举参数，错误配置尽早给出明确提示。
    if radius_style not in {"cylindrical", "tapered"}:
        raise ValueError("radius_style must be 'cylindrical' or 'tapered'.")
    # cap_style 决定圆柱端面是否封口。
    if cap_style not in {"none", "terminal", "all"}:
        raise ValueError("cap_style must be 'none', 'terminal', or 'all'.")
    # 圆周至少 6 边，否则圆柱看起来太粗糙且三角化不稳定。
    if sides < 6:
        raise ValueError("sides must be at least 6.")
    # 半径缩放必须为正，负数没有几何意义。
    if radius_scale <= 0.0:
        raise ValueError("radius_scale must be positive.")

    # 统计每个节点有多少子节点，用于判断某个 child 是否是终端节点。
    child_counts: dict[int, int] = {node_id: 0 for node_id in tree.nodes}
    # 每条父子边都让父节点的子节点计数加一。
    for parent_id, _child_id in tree.edges:
        child_counts[parent_id] += 1

    # 预先计算圆周采样角度，避免在每条边里重复生成。
    angles = np.linspace(0.0, 2.0 * np.pi, sides, endpoint=False)
    # cos/sin 数组用于把单位圆映射到每个血管段的截面平面。
    cos_values = np.cos(angles)
    sin_values = np.sin(angles)

    # points 存所有顶点；faces 存 PyVista 面片编码。
    points: list[np.ndarray] = []
    faces: list[int] = []
    # 每个面片附带一个半径和边编号，便于后续着色或调试。
    cell_radii: list[float] = []
    cell_edge_ids: list[int] = []
    # 记录退化边数量，最后在日志中提醒用户。
    skipped_edges = 0

    # 遍历每条 SWC 父子边；一条边会被转换成一个圆柱或截锥。
    for edge_id, (parent_id, child_id) in enumerate(tree.edges):
        # 取出父节点和子节点的完整记录。
        parent = tree.nodes[parent_id]
        child = tree.nodes[child_id]
        # start/end 是圆柱轴线的两个端点。
        start = parent.xyz
        end = child.xyz
        # direction 和 length 描述血管段方向与长度。
        direction = end - start
        length = float(np.linalg.norm(direction))
        # 长度过小会导致方向归一化失败，因此跳过这类退化边。
        if length <= 1e-9:
            skipped_edges += 1
            continue

        # tapered 使用父子两端各自半径，显示为截锥。
        if radius_style == "tapered":
            start_radius_um = parent.radius_um
            end_radius_um = child.radius_um
        # cylindrical 使用子节点半径表示整段半径，更贴近生成器内部边半径定义。
        else:
            start_radius_um = child.radius_um
            end_radius_um = child.radius_um

        start_radius = start_radius_um * radius_scale
        end_radius = end_radius_um * radius_scale

        # 构建垂直于血管方向的两个截面基向量。
        u, v = _basis_from_direction(direction)
        # start_ring/end_ring 保存两端圆环上每个顶点在 points 中的索引。
        start_ring: list[int] = []
        end_ring: list[int] = []

        # 沿圆周采样，分别在起点和终点生成一圈顶点。
        for c, s in zip(cos_values, sin_values):
            # normal 是当前圆周方向，位于血管截面平面内。
            normal = c * u + s * v
            # 记录起点圆环顶点索引并追加坐标。
            start_ring.append(len(points))
            points.append(start + start_radius * normal)
            # 记录终点圆环顶点索引并追加坐标。
            end_ring.append(len(points))
            points.append(end + end_radius * normal)

        # 把相邻的起点/终点圆环顶点连成四边形侧壁。
        for i in range(sides):
            # a0/a1 是起点圆环相邻两点，b1/b0 是终点圆环对应两点。
            a0 = start_ring[i]
            a1 = start_ring[(i + 1) % sides]
            b1 = end_ring[(i + 1) % sides]
            b0 = end_ring[i]
            # PyVista 四边形编码格式为 [4, p0, p1, p2, p3]。
            faces.extend([4, a0, a1, b1, b0])
            # 侧壁面片半径取两端平均值，适合 taper 模式着色。
            cell_radii.append(0.5 * (start_radius_um + end_radius_um))
            # 记录当前面片来自哪条 SWC 边。
            cell_edge_ids.append(edge_id)

        # 根节点没有父边；terminal 模式下可给入口端封口。
        add_start_cap = cap_style == "all" or (
            cap_style == "terminal" and parent.parent_id < 0
        )
        # 没有子节点的 child 是终端；terminal 模式下可给末端封口。
        add_end_cap = cap_style == "all" or (
            cap_style == "terminal" and child_counts[child_id] == 0
        )

        # 按需给起点端面加中心点和三角形扇面。
        if add_start_cap:
            center_index = len(points)
            points.append(start.copy())
            _append_cap(
                faces,
                cell_radii,
                cell_edge_ids,
                center_index,
                start_ring,
                start_radius_um,
                edge_id,
                reverse=True,
            )

        # 按需给终点端面加中心点和三角形扇面。
        if add_end_cap:
            center_index = len(points)
            points.append(end.copy())
            _append_cap(
                faces,
                cell_radii,
                cell_edge_ids,
                center_index,
                end_ring,
                end_radius_um,
                edge_id,
                reverse=False,
            )

    # 如果至少生成了一个点，vstack 得到标准 N x 3 顶点数组。
    if points:
        point_array = np.vstack(points).astype(float)
    # 如果没有可绘制边，返回空顶点数组，调用方会决定如何报错。
    else:
        point_array = np.empty((0, 3), dtype=float)

    # 将 Python 列表转换为 NumPy 数组，交给 PyVista 高效构建 PolyData。
    return MeshData(
        points=point_array,
        faces=np.array(faces, dtype=np.int64),
        cell_radii_um=np.array(cell_radii, dtype=float),
        cell_edge_ids=np.array(cell_edge_ids, dtype=np.int64),
        edge_count=len(tree.edges),
        skipped_edges=skipped_edges,
    )
