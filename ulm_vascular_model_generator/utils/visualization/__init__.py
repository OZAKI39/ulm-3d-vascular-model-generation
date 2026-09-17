"""ULM 血管 SWC 文件的可视化工具子包。"""

# MeshData/build_cylinder_mesh 负责把 SWC 树转换成圆柱 mesh 缓冲区。
from .cylinder_mesh import MeshData, build_cylinder_mesh
# planar_publication 负责平面树的论文级定量组合图。
from .planar_publication import (
    is_xz_planar,
    render_planar_publication_figure,
)
# render_swc/swc_to_polydata 负责 PyVista 渲染和 PolyData 转换。
from .pyvista_viewer import render_swc, swc_to_polydata
# swc_loader 负责读取 SWC 文本并建立节点/边结构。
from .swc_loader import SWCNodeRecord, SWCTree, list_swc_files, read_swc
# VisualizationConfig/load_visualization_config 负责解析可视化 YAML。
from .yaml_config import VisualizationConfig, load_visualization_config

# __all__ 声明可视化子包推荐对外使用的 API。
__all__ = [
    # mesh 构建结果和函数。
    "MeshData",
    # SWC 解析结果类型。
    "SWCNodeRecord",
    "SWCTree",
    # 可视化配置对象。
    "VisualizationConfig",
    # 可视化和文件读取函数。
    "build_cylinder_mesh",
    "is_xz_planar",
    "list_swc_files",
    "load_visualization_config",
    "read_swc",
    "render_swc",
    "render_planar_publication_figure",
    "swc_to_polydata",
]
