# 源码复用审计

| 官方项目 | 实际 Git commit | 许可 |
|---|---|---|
| [OpenPIV](https://github.com/OpenPIV/openpiv-python) | 91e88a8ee9289d6c4ca79722d9759c5258abd31d | GPL-3.0 |
| [TorchPIV](https://github.com/NikNazarov/TorchPIV) | 17ae4c1a4d4950efecfdfd9c404d77036a6f8c80 | MIT |
| [CuPy PIV](https://github.com/ali-sh-96/openpiv-python-gpu) | d50d505109f9ef9f8ddbd971380ee2ac74f38682 | GPL-3.0 |

FFT 相关、峰值定位、SNR、异常向量验证、多遍形变与 GPU 批量内核均直接调用上游；没有另写 FFT 或峰值求解器。OpenCV 承担候选/积分图/配准，SciPy 承担原六参数鲁棒拟合及矩阵对数，PyAV/FFmpeg 复用 F10V 编码和完整解码检查。GPL 组件不是可重新标为 MIT 的代码；分发时必须保留原许可与对应源码要求。许可证全文在 evidence/licenses。

实际兼容改动：TorchPIV 包导入强制 GUI，薄 namespace 适配绕过 GUI 初始化后调用未改动 PIVbackend；OpenPIV local_median_val 并不支持 backend 参数，按真实 API 调用；CuPy cached s2n_ratio 属性返回 None，读取其已计算 corr.sig2noise；CuPy bool mask cumsum 默认为 int64，而 kernel 读 int32，改为 `cp.cumsum(mask, dtype=self.dtype_i)`。最后一项属于 ABI 修复，不修改相关算法，diff/哈希保存在 evidence/cupy_mask_compat。修复后 B2/B3/B4 三种选窗掩码在共同窗口上的位移差均为 0。

OpenPIV 0.26.1 pip wheel 本身未包含可用 Rust 扩展，使用官方 maturin 构建后显式导入 openpiv_rust 并执行 backend='rust'；没有自动回退。官方 CuPy tutorial test0 实际执行，12600 个有限向量、20 个替换/拒判仅作为示例，不混入正式观测。

FluidImage 和 CPU FFTW/Numba 备用方案 NOT_RUN：Rust 已正常工作，GPU PIV 核心不再是主要瓶颈，本轮不迁移更多框架。配置与版本锁见 evidence/remote_setup 和 evidence/local_requirements_final.lock。
