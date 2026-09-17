"""
SWC 读取工具。

生成器输出的是文本 SWC 文件；可视化前需要把每一行解析成节点，并根据
parent_id 建立边列表。这个模块只负责读文件和基础校验，不负责画图。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SWCNodeRecord:
    """解析后的一行 SWC 节点记录。"""

    # SWC 第一列：节点编号，同一个文件内必须唯一。
    node_id: int
    # SWC 第二列：节点类型；本项目主要用于兼容标准 SWC 格式。
    node_type: int
    # SWC 第三到第五列：节点三维坐标，单位沿用生成器输出的微米。
    xyz: np.ndarray
    # SWC 第六列：节点半径，用于可视化时生成圆柱截面。
    radius_um: float
    # SWC 第七列：父节点编号；负数表示根节点没有父节点。
    parent_id: int


@dataclass(frozen=True)
class SWCTree:
    """SWC 节点字典和父子边列表。"""

    # 原始 SWC 路径，便于错误信息和渲染标题显示来源。
    path: Path
    # 用节点编号作为 key，便于 O(1) 查找父节点或子节点。
    nodes: dict[int, SWCNodeRecord]
    # 每条边保存为 (parent_id, child_id)，后续会把每条边画成一个圆柱。
    edges: list[tuple[int, int]]

    @property
    def node_count(self) -> int:
        """返回节点总数，供日志和测试使用。"""
        return len(self.nodes)

    @property
    def edge_count(self) -> int:
        """返回父子边总数；一条边对应一个可绘制的血管段。"""
        return len(self.edges)

    @property
    def coordinates(self) -> np.ndarray:
        """把所有节点坐标堆叠成 N x 3 数组，方便做包围盒或统计。"""
        # 空树没有坐标，返回形状正确的空数组，避免调用方特殊处理 None。
        if not self.nodes:
            return np.empty((0, 3), dtype=float)
        # vstack 将若干个长度为 3 的坐标向量拼成二维数组。
        return np.vstack([node.xyz for node in self.nodes.values()])

    @property
    def radii_um(self) -> np.ndarray:
        """把所有节点半径取出成一维数组，方便做半径着色或统计。"""
        # 空树同样返回空数组，保持 API 类型稳定。
        if not self.nodes:
            return np.empty((0,), dtype=float)
        # 明确 dtype=float，避免整数半径输入导致后续数值运算类型不一致。
        return np.array([node.radius_um for node in self.nodes.values()], dtype=float)


def list_swc_files(directory: Path) -> list[Path]:
    """列出目录中的 SWC 文件。"""
    # Path.glob 只查当前目录；按小写文件名排序可让输出顺序稳定。
    return sorted(Path(directory).glob("*.swc"), key=lambda path: path.name.lower())


def read_swc(path: Path) -> SWCTree:
    """读取标准 7 列 SWC：id type x y z radius parent。"""
    # 统一转为 Path，允许调用方传入字符串或 Path。
    path = Path(path)
    # nodes 用字典存储，后面检查 parent_id 是否存在时可快速查找。
    nodes: dict[int, SWCNodeRecord] = {}
    # edges 暂时为空；等所有节点读完后再建立边，避免父节点出现在子节点之后时报错。
    edges: list[tuple[int, int]] = []

    # SWC 由本项目以 ASCII 写出；按 ASCII 读取可尽早暴露异常字符。
    with path.open("r", encoding="ascii") as handle:
        # line_number 用于在报错时指出具体哪一行不合法。
        for line_number, line in enumerate(handle, start=1):
            # 去掉首尾空白，方便判断空行和注释行。
            stripped = line.strip()
            # 空行或 # 开头的注释行不含节点数据，直接跳过。
            if not stripped or stripped.startswith("#"):
                continue
            # 标准 SWC 以空白分隔字段。
            fields = stripped.split()
            # 本读取器至少需要前 7 列，少列说明文件不是合法 SWC。
            if len(fields) < 7:
                raise ValueError(f"{path}:{line_number}: expected 7 SWC columns.")

            # 前两列是整数编号和类型。
            node_id = int(fields[0])
            node_type = int(fields[1])
            # 坐标列转成浮点数，保持亚像素/微米精度。
            x, y, z = (float(fields[2]), float(fields[3]), float(fields[4]))
            # 半径必须为正，否则圆柱 mesh 无法构造。
            radius_um = float(fields[5])
            # 父节点编号决定血管段连接关系。
            parent_id = int(fields[6])

            # 节点编号重复会导致父子查找歧义，因此直接报错。
            if node_id in nodes:
                raise ValueError(f"{path}:{line_number}: duplicate node id {node_id}.")
            # 半径为 0 或负数时，圆柱退化或方向不明确，提前阻止。
            if radius_um <= 0.0:
                raise ValueError(f"{path}:{line_number}: radius must be positive.")

            # 将这一行整理成不可变记录，避免可视化过程中误改原始数据。
            nodes[node_id] = SWCNodeRecord(
                node_id=node_id,
                node_type=node_type,
                xyz=np.array([x, y, z], dtype=float),
                radius_um=radius_um,
                parent_id=parent_id,
            )

    # 所有节点已读入后，再把 parent_id 转成边列表。
    for node in nodes.values():
        # parent_id < 0 表示根节点，不生成父子边。
        if node.parent_id < 0:
            continue
        # 如果父节点不存在，说明 SWC 文件内部引用断裂。
        if node.parent_id not in nodes:
            raise ValueError(f"{path}: node {node.node_id} references missing parent {node.parent_id}.")
        # 保存父子关系；可视化时一条边会变成一个圆柱。
        edges.append((node.parent_id, node.node_id))

    # 按子节点编号排序，让 mesh 构建和调试输出更稳定。
    edges.sort(key=lambda edge: edge[1])
    # 返回统一的树对象，调用方不需要直接处理原始文本行。
    return SWCTree(path=path, nodes=nodes, edges=edges)
