# 新增和修改文件清单

生产代码版本：`98d705f13babdfa0cba6c12331623292a588cb77`。完整逐文件路径、大小及 SHA256 见 `FILE_INVENTORY.csv`。

原先117个物理模块字节保持不变。本轮只新增 population_v1 包、专用测试、BraVa population_v1 目录；不改 CFD、旧源分布、旧轨迹/动画。开发过程修复限定在新增模块。

|模块|用途|
|---|---|
|`environment.py`|权威环境与哈希/端口适配|
|`source.py`|精确流量源、事件 RNG、单次准入|
|`state.py`|CPU SoA、稳定ID|
|`engine.py`|CPU 共享时钟与群体守恒|
|`reference.py`|原物理函数的私有依赖绑定|
|`pairs.py`|精确球对间隙/候选/新交互契约|
|`components.py`|局部组件与稀疏耦合接触|
|`safety.py`|继承数值保护和原支持证据|
|`output.py`|完整接受时间状态和账本|
|`run.py`|CLI、真实覆盖停源、删失|
|`gpu_engine.py`|实际GPU动力学及明确CPU小子集回退|
|`gpu_field.py`|静态FEM上传、tetra walking|
|`gpu_wall.py`|真实有限壁面GPU查询/路径|
|`gpu_neighbors.py`|GPU空间散列及精确gap|
|`gpu_resistance.py`|GPU孤立泡阻力和接触|
|`gpu_state.py`|GPU SoA|
|`gpu_safety.py`|GPU进展历史和终态CPU证据|
|`gpu_batch.py`|各自时钟的多组数值叠批|
|`gpu_extension.py`|FP64 CUDA编译|
|`bvh.py`|静态BVH和tetra邻接|
|`provenance.py`|规范JSON和内容哈希|
|`cuda/kernels.cu`, `cuda/bindings.cpp`|FEM/散列/有限三角形 CUDA kernels 与 PyTorch 绑定|

测试分 Phase2/3/4/6/7/8 和 safety/output；Phase5 真实点集及实际场验证为 scripts 下独立执行脚本。绘图、报告、校验脚本均在 `scripts/`；未将上述工具散落到原正式目录。

已终止的开发失败输出保持原路径和原字节。实际环境二进制、CUDA编译缓存和大体网格属于服务器环境/原项目，交付清单提供路径与哈希，不复制到新增源码包。
