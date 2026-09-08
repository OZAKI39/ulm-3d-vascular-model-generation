# 本轮复用与来源

所有来源项目只读。本轮没有从旧工程 `sys.path` 导入或启动完整 CFD pipeline。精确源文件 SHA256、Mirheo 提交与编译库身份见正式结果包 `provenance.json`；整套已有文件的保护快照位于 `test_code/outputs/fluid_physics/setup_20260908T154306_330273Z/protected_before.json` 及交付后的完整性审计。

|能力|已有来源|复用方式|
|---|---|---|
|血管身份、数组、端口、单位|`py_scripts/vessel_geometry/export.py::load_package`|直接调用已验收加载器，验证所有包内哈希，不重建形状|
|路径防护、JSON、SHA256|`py_scripts/vessel_geometry/io.py`|直接复用 `safe_output`、`require_file`、`read_json`、`write_json`、`sha256_file`|
|旧纯流体工况|旧工程 `configs/cfd_flow.yaml`、接受 run 的生成 Lua / final report / runtime contract|独立读取，校对常量、rho Q、压力偏置和源身份|
|旧物理截面|`physical_port_flux_plane_contract_v3.json`、旧 `physical_port_flux.py`、冻结刚体变换|保留有限孔口 UV 轮廓；仅逆变换位置和方向元数据|
|旧验收|旧 `utils/cfd_flow/validated_contract.py`、`steady_state.py`、接受 run 的 `reference_scaled_base_physical_flux_history.json`|重新实现少量数学纯函数；从实际数据验证分流与密度流量分母；不移植求解器 I/O 或 LBM 检查|
|GPU 环境与 MPI|现有 `scripts/activate_mirheo.sh`、`scripts/run_official_case.sh`|source 原激活脚本；原 OpenMPI `--bind-to none -np 2`，一个计算 rank、一个后处理 rank|
|小周期液体施力|`vendor/Mirheo/tests/flow/double_poiseuille.py`、`core/integrators/forcing_terms/periodic_poiseuille.h`|改编接口调用；核对每粒子力、x 方向、y 半盒变号和质量因子|
|DPD 随机耗散对应|`core/interactions/pairwise/kernels/dpd.h`|只读核实 sigma、wr、dt 与质量，不重复添加随机力|
|压力|`tests/stress/pressure.py`、`plugins/virial_pressure.cu`、`kernels/stress_wrapper.h`|改编原生插件接线；补上动能和体积，不重复配对因子|
|温度、时间标签|`plugins/stats.cu`、`core/simulation.cpp`|读取真实 CSV，核对均动量含义和同回调两种时间标签；COM/分箱热速度统计由本轮实现|
|入口出口设计|`bindings/plugins.cpp`、`plugins/velocity_inlet.cu`、`density_control.cu`、`outlet.cu`|记录准确原生能力和限制，未输出虚构 setter 或压力储液 API|
|3D/HTML|`test_code/review_vessel_geometry.py::make_figure`、`open_windows_file`|复用完整几何视图和 Windows 打开方式；新增真实统计曲线和冻结截面线|
|浏览器检查|`test_code/check_vessel_geometry_browser.cjs`|改编为新的独立 CDP checker，保留原文件；无 npm 依赖和 Web 服务|

Mirheo 版本为本地提交 `8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf`。示例接口改编遵循原项目 MIT 许可证（Copyright 2019 ETH Zurich）；每个 GPU 任务均保存原样 `MIRHEO_LICENSE.txt`，上游源码已有 Copyright 注释全部保留。没有改写或重编译 vendor。

旧用户工程根目录没有发现一份可据以宣称整个项目采用同一许可证的 LICENSE。这里仅依据用户授权读取本地材料并重新实现数学定义，不把第三方子目录许可证错误套到整个旧工程，也不新增覆盖它的许可证声明。

来源字段类别使用 `SOURCE_VALUE / USER_ASSUMPTION / DESIGN_CHOICE / MEASURED / DERIVED / PROPOSED / UNVERIFIED`。准备包的早期诊断映射曾补充为带接受流量采样证据的新不可覆盖准备包，源物理值、GPU 输入与预算没有变化，父版本保留在 `prepared_package.json` 的 `supplements` 链中。
