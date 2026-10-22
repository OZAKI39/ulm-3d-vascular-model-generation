# s1-8 显式血管支撑的开源来源审计

当前正式方向为 **CUSTOM_LITERATURE_AWARE_SUPPORT**。血管、半径、分叉、端口、打印姿态和 2 mm 五面盒全部冻结。原 Bambu、Orca、Guided Orca 结果保留为 **BENCHMARK_ONLY**。Bambu 在正式验证中关闭自动支撑，仅对已经生成的显式支撑和生产模型切片。

## FreeFloating：保留了什么

实际检查并固定 [regrow1123/FreeFloating](https://github.com/regrow1123/FreeFloating) 到 `5462f823abf01661e9abeb85d9bc871e01d4ec39`。该版本为 MIT，原作者 Seongje Jang，2019。`support_freefloating.py` 移植相邻射线入口/出口有向图、层号量化及加权图覆盖，使用 trimesh/Embree 替代 OpenGL 深度链表。原许可证、版权和检查过的源文件保留在 `third_party/support_research/FreeFloating/`。

原始应用的编译探针在 Boost 头文件处失败，日志保留；没有宣称运行过完整原始程序。永久测试验证图规则、量化、局部低点和真实射线。生产检测保留 `tables/freefloating_raw_points.csv`，再把有效证据映射到命名分支和弧长，不能将每个原始点直接变成打印支撑。

## Clever Support：作者下载包实际包含什么

实际访问了 [作者实验室论文页](https://cgvlab.github.io/cgvlab/www/publications/Vanek14SGP/)。页面中的大写 `.ZIP` 链接返回 404；检查 [作者服务器目录](https://www.cs.purdue.edu/homes/bbenes/papers/) 后，从小写后缀下载了 [Vanek14SGP.zip](https://www.cs.purdue.edu/homes/bbenes/papers/Vanek14SGP.zip)。

文件为 6,559,878 字节，SHA256 为 `91d9c031bf0e2b5a37c444debaa20d603033b099c8413d5cfa3a444901e29afe`。实查 ZIP 只有 14 个 OBJ：7 个原始模型和 7 个带支撑模型；**没有程序源码、构建文件或代码许可**。因此没有可编译的“官方代码应用”，也没有把这些模型当作 MIT 源码复制。下载包、目录清单和论文留在本地外部审计缓存，项目保存获取记录与摘要。

实际阅读论文第 5、6 节：从高到低处理候选点、合并向下支撑锥体、建立支撑图，并依据长度和角度考虑支杆稳定性。当前适配与论文有明确区别：不重新搜索打印姿态，不允许以血管或侧壁为根，仅允许盒内底板的受控底脚；不把论文中的原打印机经验参数当成 X1C ABS 的强度规律。

## meshy：真实可审计的 MIT 实现参考

[meshy](https://github.com/0x00019913/meshy) 固定到 `4767c520d99c7c6275b1cb0d8fd843e548eaac53`，MIT，Anton Bovin，2017。实际读取 `js/supportGenerator.js` 的高点优先处理、锥体配对、碰撞回退、累积载荷和渐细支杆，另外检查 `js/calculate.js`。保留原 LICENSE 和源文件。

当前 Python 模块沿用其 `r + k sqrt(load)` 加粗规律，并据所审计的树合并流程进行显式几何适配。几何实体由 trimesh/manifold 生成，碰撞和真实三角形间隙交给 FCL；没有复制浏览器 UI、WebGL、G-code 生成器。新增了完整尖端间隙校正、可折断窄颈、只在盒内底板落根、自身逐层支撑检查，以及折颈后整组件的扫掠移除检查。

FreeFloating 仓库中的 `sources/Vanek/` 实际是顶点、边和面悬垂检测的比较实现，并非完整 Clever 树生成器；审计明确区分这两者，没有把它当作原作者完整实现。

## 外部切片器与短桥论文

[OrcaSlicer](https://github.com/OrcaSlicer/OrcaSlicer) 审计源码版本为 `c023632a7d7204a163dbcba7a2972581034452b3`，AGPL-3.0；实际外部程序为官方 2.4.2，发布摘要和本地运行库摘要已保存。[CuraEngine](https://github.com/Ultimaker/CuraEngine) 审计版本为 `741be28fdfecd4428b84ab6046098e611e10dd60`，AGPL-3.0，仅阅读树支撑相关实现，没有用于本次实际生产切片。两者源码留在外部审计缓存，没有大段复制进项目。

Orca 原生树和 Guided Orca 均只作为基准。Bambu 原生树也只作基准，但 Bambu X1 Carbon 0.2 mm 的实际本地程序承担显式装配、独立支撑及 coupon 的 **Auto Support OFF** 验证。工艺值从本机保存的 `0.10mm Standard Vessel` / ABS Vessel 解析，不生成自己的挤出策略。历史 Guided Orca 的 3MF schema 已经实际引擎读取与切片验证，原结果仍完整保留。

[Bridging the Gap 作者页面](https://www.jdumas.org/publication/2014/scaffoldings/) 提供立柱加短桥的算法参考，并提及 IceSL 软件；本项目没有取得并调用可嵌入的官方支撑源码。独立短桥候选只标记为 `GEOMETRY_CANDIDATE`；4/6/8 mm 的研究结果保留，不能替代真实桥接试片。生产修复优先采用加粗、合并或竖直立柱，不能用新增血管接触掩盖不稳定。

所有来源、实际检查文件、复用内容和改动见 `third_party/support_sources.json`。软件仅执行本地切片、预览与导出，不上传或启动打印。软件输出不代表 ABS 强度、颈部断裂力、真实支撑稳定性或拆除成功已经验证。
