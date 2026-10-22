# 官方来源、权重与复用

EXPLORATORY · PHYSICS_UNCALIBRATED · training_ready=false

## LRAFT-ST

官方仓库：https://github.com/wuwuwuas/LRAFT-ST.git

固定提交：`5ebe4bf9e769ebd4a26a6d2284753cc79b2cbb10`；许可：`NO_LICENSE_FOUND_DO_NOT_REDISTRIBUTE`。

## co-tracker

官方仓库：https://github.com/facebookresearch/co-tracker.git

固定提交：`82e02e8029753ad4ef13cf06be7f4fc5facdda4d`；许可：`FOUND_READ_ATTACHED`。

## UnLiteFlowNet-PIV

官方仓库：https://github.com/erizmr/UnLiteFlowNet-PIV.git

固定提交：`47cf7d3423aa48b69ee49a86067176a0f68a4174`；许可：`NO_LICENSE_FOUND_DO_NOT_REDISTRIBUTE`。

CoTracker LICENSE.md为CC-BY-NC 4.0；本地科研验证，禁止据此假定商用许可。LRAFT和UnLite未发现LICENSE，授权不明；第三方源码/权重没有纳入本项目Git交付。POT为MIT，直接调用官方求解器。OpenPIV、pycpd、F10U几何转换器均复用。

- `/workspace/f10w_vision_20261017/models/lraft_student_dataset1.tar`：11708931 bytes，SHA-256 `2ce5026afee3a8a0a535a472591fa73abcf119efd07cb07da45b4809aa12c36f`。
- `/workspace/f10w_vision_20261017/models/cotracker3_scaled_online.pth`：101695610 bytes，SHA-256 `205d34789f19699d64b22cf93f9b697f15f28d4025240e31532e504109837218`。
- `/workspace/f10w_vision_20261017/third_party/UnLiteFlowNet-PIV/models/UnsupervisedLiteFlowNet_pretrained.pt`：86203349 bytes，SHA-256 `96ab5270dd27c4cc5002b26cc8a1ed45d72df1114b38c4bdeee802e806c0bb16`。

LRAFT权重来自官方README链接的Google Drive学生Dataset1（ID 18tTvU8WjquJLY84GHrxrgeHQXpfZ5UQh）；CoTracker来自官方Hugging Face facebook/cotracker3/scaled_online.pth；UnLite使用官方仓库自带权重。下载脚本和日志保留。UnLite仅将CuPy旧memoize/compile_with_cache接口适配为memoize/RawModule(backend=nvcc)，CUDA核内容、网络和权重未改变；原失败、diff、before/after SHA保存在evidence/unlite_compat。
