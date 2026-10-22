# 活动窗口与批处理

每帧全通道快速候选+纹理+差分，不依赖旧目标附近。安全且有候选的静止区域仍激活，新出现窗口立即更新；高候选密度>=2也每帧。低活动细 PIV比较1/3/6帧周期，无效/跳过速度为null。

CPU Rust和TorchPIV用官方 sliding_window_array 批量gather成无重叠条带，再调用原完整求解器；CuPy使用官方field_mask，经ABI兼容修复后工作。B2/B3/B4 all-safe/activity/checker三种mask与规则计算共同窗口位移最大差0。

真实部分审核505个中心：466/505在8px内有检测，503/505被活动窗覆盖。初版几何漏掉的真实泡已恢复。合成新位置 probe：静止、小位移、稀疏及dropout中心均被覆盖；boundary漏1个；低对比只覆盖10/120。静态PDMS假候选会出现，几何排除与后续PIV gate控制，但无法证明剩余所有候选都是泡。真实PZT-OFF状态无元数据，NOT_IDENTIFIED。

活动筛选验证仅PARTIAL，不设全通道召回率已验证。规则GPU批处理在相同图像pilot比稀疏gather/mask更快，因此采用规则主线；窗口减少不等于完整管线更快。数据在 channel_audit、post_freeze_diagnostics及channel_pilot_v2。
