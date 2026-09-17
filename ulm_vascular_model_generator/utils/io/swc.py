"""
把内部 Vessel 树导出为 SWC 格式。

内部算法用 `Vessel` 表示“边”（一段血管圆柱），而 SWC 格式用“节点”
表示树：每一行是一个点，最后一列指向父节点。为了兼容常见神经/血管
可视化工具和下游仿真，这里把每条 vessel 的远端点写成一个 SWC 节点。
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core.models import Vessel


@dataclass
class SWCNode:
    """SWC 文件中的一行节点。"""

    # SWC 第一列：节点编号。
    node_id: int
    # SWC 第二列：节点类型，常见工具会用它区分 soma/axon/dendrite/vascular 段。
    node_type: int
    # SWC 第三到第五列：节点坐标。
    xyz: np.ndarray
    # SWC 第六列：节点半径。
    radius_um: float
    # SWC 第七列：父节点编号；-1 表示无父节点。
    parent_id: int


def vessels_to_swc(vessels: list[Vessel],
                    root_type: int = 2,
                    other_type: int = 3) -> list[SWCNode]:
    """
    把 Vessel 列表转换成 SWCNode 列表。

    root vessel 会写两个节点：近端 root 节点和远端节点。
    之后每条子 vessel 只写它的远端节点，父节点指向父 vessel 的远端 SWC
    节点。这样 SWC 拓扑就和 Vessel 树拓扑一致。
    """
    # 空 Vessel 树对应空 SWC。
    if not vessels:
        return []

    # parent vessel id -> child vessel ids。这里复制一份，避免后续操作改动原树。
    # children_map 用于广度优先遍历 vessel 树。
    children_map: dict[int, list[int]] = {v.vid: list(v.children) for v in vessels}
    # vessel_distal_swc 记录每条 vessel 的远端点对应哪个 SWC 节点。
    vessel_distal_swc: dict[int, int] = {}
    # swc 保存最终输出的节点列表。
    swc: list[SWCNode] = []
    # root vessel 是 parent_id < 0 的 vessel。
    root_ids = sorted(v.vid for v in vessels if v.parent_id < 0)
    # 防御性处理：如果没有 root，就把第一条 vessel 当作 root。
    if not root_ids:
        root_ids = [vessels[0].vid]

    # 每个 root 独立写一棵 SWC 子树。
    for root_id in root_ids:
        root = vessels[root_id]
        # root vessel 的近端点需要单独写成 SWC 根节点。
        proximal_id = len(swc)
        swc.append(SWCNode(
            node_id=proximal_id,
            node_type=root_type,
            xyz=root.x_p.copy(),
            radius_um=root.radius,
            parent_id=-1,
        ))

        # root vessel 的远端点写成第二个节点，父节点指向近端 root 节点。
        distal_id = len(swc)
        swc.append(SWCNode(
            node_id=distal_id,
            node_type=other_type,
            xyz=root.x_d.copy(),
            radius_um=root.radius,
            parent_id=proximal_id,
        ))

        # 保存 root vessel 远端对应的 SWC 节点 id，子 vessel 会把它作为父节点。
        vessel_distal_swc[root.vid] = distal_id

        # 从 root 的孩子开始广度优先遍历整棵 vessel 子树。
        queue: deque[int] = deque(children_map[root.vid])
        while queue:
            # 取出下一条 vessel。
            vid = queue.popleft()
            v = vessels[vid]
            # 当前 vessel 的 SWC 父节点，是父 vessel 的远端 SWC 节点。
            parent_swc = vessel_distal_swc[v.parent_id]
            # 当前 vessel 只需要写远端点，近端点已经由父 vessel 的远端表达。
            node_id = len(swc)
            swc.append(SWCNode(
                node_id=node_id,
                node_type=other_type,
                xyz=v.x_d.copy(),
                radius_um=v.radius,
                parent_id=parent_swc,
            ))
            # 保存当前 vessel 远端节点 id，供它的孩子引用。
            vessel_distal_swc[vid] = node_id
            # 将当前 vessel 的孩子加入队列继续遍历。
            queue.extend(children_map[vid])

    # 返回 SWC 节点列表，写文件函数负责具体文本格式。
    return swc


def write_swc(nodes: list[SWCNode], path: Path) -> None:
    """按标准 SWC 文本格式写文件。"""
    # 确保输出目录存在。
    path.parent.mkdir(parents=True, exist_ok=True)
    # SWC 只包含 ASCII 数字和空格；newline="\n" 保证跨平台输出稳定。
    with path.open("w", encoding="ascii", newline="\n") as handle:
        # 每个 SWCNode 写成一行 7 列。
        for node in nodes:
            # 拆出坐标，便于格式化。
            x, y, z = node.xyz
            # 半径和坐标保留 6 位小数，兼顾精度和文件体积。
            handle.write(
                f"{node.node_id:d} {node.node_type:d} "
                f"{x:.6f} {y:.6f} {z:.6f} "
                f"{node.radius_um:.6f} {node.parent_id:d}\n"
            )


# ======================================================================
# Summary
# ======================================================================


def summarize(nodes: list[SWCNode]) -> str:
    """生成一行简短统计信息，供 CLI 打印。"""
    # 空树没有几何统计，直接返回固定文本。
    if not nodes:
        return "(empty tree)"
    # 堆叠所有坐标，用于计算空间范围。
    coords = np.vstack([n.xyz for n in nodes])
    # 收集所有半径，用于输出 min/median/max。
    radii = np.array([n.radius_um for n in nodes])
    # parent_ids 用来判断哪些节点被其他节点作为父节点。
    parent_ids = np.array([n.parent_id for n in nodes])
    # has_child 是至少有一个孩子的 SWC 节点集合。
    has_child = set(parent_ids[parent_ids >= 0].tolist())
    # 没有孩子的节点就是叶子节点。
    leaves = sum(1 for n in nodes if n.node_id not in has_child)
    # extent 是 x/y/z 三个方向的包围盒尺寸。
    extent = coords.max(axis=0) - coords.min(axis=0)
    # 返回紧凑摘要，便于命令行快速确认生成结果。
    geometry_text = (
        f"plane_y_um={coords[0, 1]:.1f}"
        if extent[1] <= 1.0e-6
        else "geometry=volumetric_3d"
    )
    return (
        f"nodes={len(nodes)}, leaves={leaves}, "
        f"extent_um=({extent[0]:.1f}, {extent[1]:.1f}, {extent[2]:.1f}), "
        f"{geometry_text}, "
        f"radius_um=({radii.min():.2f}, {float(np.median(radii)):.2f}, {radii.max():.2f})"
    )
