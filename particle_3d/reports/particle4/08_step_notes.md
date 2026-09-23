![08_pair_physical_time_subdivision](figures/08_pair_physical_time_subdivision.png)

应该看什么：看端点分离的大步是否仍检测到中途相撞。

实际看到什么：三个请求步长均在 0.375 秒接触，并完整处理到 1 秒。

有没有异常：所有细分区间都有真实起止时间，失败尝试不消耗时间。
