# VascularMD 第三方依赖记录

- 来源：<https://github.com/megdec/vascularmd>
- 固定提交：`770feb8fbfb591d6d43f966db57b1465010b9a13`
- 位置：`third_party/vascularmd/`，保留官方源码和 GPL-3.0 LICENSE。
- 本地补丁：[vascularmd_compat.patch](vascularmd_compat.patch)。

仅有以下兼容修改：

1. `ArterialTree.py`、`Spline.py`、`Model.py`、`Nfurcation.py`、`utils.py` 的内部导入改为相对导入，并添加 `__init__.py`。官方使用 `from utils import *`，会与本项目同名 `utils` 包发生冲突；包隔离避免污染 `sys.modules` 或覆盖已有工具。
2. `Nfurcation.smooth()` 的 `mesh.smooth(n_iter, ...)` 改为 `mesh.smooth(n_iter=n_iter, ...)`，适配 PyVista 的关键字参数要求。迭代次数、松弛系数、算法和原生回投影流程不变。

没有改动 AIC、控制点求解、独立空间/半径拟合、分叉重建、合并分叉或网格数学模型。论文为 Decroocq et al., *Medical Image Analysis*, 89 (2023), 102912，DOI: [10.1016/j.media.2023.102912](https://doi.org/10.1016/j.media.2023.102912)。本次依据用户提供的本地原文核对方法及默认参数。

在没有第三方目录的全新工作区，从项目根目录复现：

```bash
git clone https://github.com/megdec/vascularmd.git third_party/vascularmd
git -C third_party/vascularmd checkout 770feb8fbfb591d6d43f966db57b1465010b9a13
git -C third_party/vascularmd apply ../vascularmd_compat.patch
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python -m pip install -r requirements-vascularmd.txt
```

已经应用补丁的工作区不要再次应用。可核查：

```bash
git -C third_party/vascularmd rev-parse HEAD
git -C third_party/vascularmd diff
git -C third_party/vascularmd apply --reverse --check ../vascularmd_compat.patch
```

每份 QC 中均记录实际提交和五个核心源文件的 SHA-256；[vascularmd.lock.json](vascularmd.lock.json) 保存本次最终集成版本摘要。
