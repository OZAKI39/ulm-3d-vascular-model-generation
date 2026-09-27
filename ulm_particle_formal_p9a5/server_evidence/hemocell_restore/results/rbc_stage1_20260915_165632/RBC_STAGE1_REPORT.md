# RBC_ONLY_VALIDATION_STAGE_1：当前生产版本

**RBC_ONLY_VALIDATION_STAGE_1 = FAIL。** 停在真实RBC网格几何适配硬门。MPI1 sanity及CASE A/B/C全部未运行，RBC_TIMESTEPS=0。另有当前Guo血管配置缺少RBC-wall作用实现，以及固定名义缩放的实际离散体积偏离50 µm³两项需审阅的问题。

## 本轮实际完成的工作

在Vast RTX4090实例建立独立工作目录/workspace/hemocell_restore/work/rbc_stage1_20260915_165632。从只读minimal restore恢复并SHA256验证1512个原生文件；完整构建HemoCell CPU/MPI静态库和新的原生RBC网格生成程序。没有使用旧服务器未知binary，没有修改HemoCell core或纯流体基线。

HemoCell commit=5a410848bd5c57d5ae1c171112e78eab4a82e650；Palabos commit=05712164d940a42e06afdd705249912fa0c49f14。额外从已校验Palabos原始源码包应用既有HemoCell patch重建，与1421个Palabos文件逐一一致。编译器、MPI、CMake、并行HDF5、flags与SHA256见RBC_STAGE1_BUILD_PROVENANCE.json。并行HDF5配置与链接通过；RBC HDF5运行输出管线未进行测试。

构建范围必须区分：原生CPU库和geometry construction executable为PASS；没有发布或声称完成可执行IBM timestep的生产case binary。几何硬门失败后，不继续建立依赖有效spawn的耦合case入口，也没有MPI1 sanity或MPI4运行证据。原生网格生成在服务器完成（0 timestep），几何搜索和独立复核在WSL已有NumPy/VTK环境完成。

## 介质合同保持

先前PBS/BSA pure-fluid smoke与SUSPENDING_MEDIUM_NUMERICS=PASS保留；ρ=1000 kg/m³、ν=1e-6 m²/s、tau=1、dt=6.6599470814617197e-09 s。三个出口rho_LU为1.0000484343922278、1.0004402376508774、0.99995437727568648。直接复用已经验证的生成合同，prepare_numerics仍是唯一生成源，本轮不再生成另一套dt或出口密度。副本与生成脚本、输入、回执的SHA256相符。

本轮没有RBC耦合solver，因此新的RBC solver实际读取dt/出口密度一致性是NOT_RUN，不能借用纯流体PASS来冒充RBC运行PASS。介质仍为1×PBS+1%BSA、25°C的开发假设，EXPERIMENTALLY_MEASURED=NO，最终实验合同PENDING。入口倍率1.1197286861799598未重新校准，新介质校准仍NOT_PERFORMED。

## 几何失败证据

细节见RBC_GEOMETRY_FIT_REPORT.md。固定scale=0.8220706914434901，名义直径6.428592807µm。原生网格642顶点/1280三角形，实际体积45.046330078320µm³（比名义50低约9.907%），面积87.307010647669µm²。该差异是初始化几何差异，不是运行volume error。

8974155个确定性粗筛姿态及预定连续细化未找到安全初态。最佳失败姿态有154个顶点在lumen外、最深越界约0.597581µm，并有186个曲面接触对。未自动缩小RBC，也未调整网格、kLink/kArea/kBend/kVolume、粘度或tau。模型常数15/5/80/20仅被冻结为本轮开发定义；膜力、IBM、体积与面积时间演化均未验证。

## CASE B 的独立阻塞

HemoCell库确实存在enableBoundaryParticles、populateBoundaryParticles、applyBoundaryRepulsionForce。现有血管实现没有启用该路径。其壁面粒子生成依赖Dynamics::isBoundary()，而当前血管是BGK/NoDynamics配合Guo off-lattice processors，不提供该格点壁面粒子集合。流体Guo no-slip不能视作RBC排斥力。

所以WALL_INTERACTION_IMPLEMENTATION=ABSENT（限定当前production vascular setup），CASE_B=BLOCKED_WALL_INTERACTION_IMPLEMENTATION。未臆造force/cutoff，未把Guo换成bounce-back，也未修改HemoCell核心来绕过此条件。库中接口存在不等于当前血管应用具备已验证的壁面相互作用。

## 阶段与未验证项

|阶段|结果|实际步数|
|---|---|---:|
|Native CPU/MPI library + reference mesh build|PASS|0|
|真实RBC几何fit|FAIL|0|
|MPI1 sanity|NOT_RUN_GEOMETRY_GATE_FAILED|0|
|CASE A / MPI4|NOT_RUN_GEOMETRY_GATE_FAILED|0|
|CASE B / MPI4|BLOCKED_WALL_INTERACTION_IMPLEMENTATION，同时几何门未通过|0|
|CASE C / MPI4|NOT_RUN_GEOMETRY_GATE_FAILED；未选定2–4个安全位置|0|

没有以MPI1回退替代MPI4。没有RBC count/ownership、IBM耦合、nonfinite时间序列、force/velocity runaway、面积/体积漂移、动态triangle inversion或运行中穿墙的PASS证据。对应数值为UNVERIFIED/NOT_RUN，不能填0冒充验证。wall metric仅静态验证；runtime wall penetration为UNVERIFIED。

RBC_STAGE1_CONTRACT.json记录用户指定阈值、计划步数/核数、监测与输出政策，以及因无安全spawn而故意保留的null位置。任何未来运行都必须先解决几何定义/适配与壁面实现问题，重新冻结可审阅的有效初态和运行合同；本轮失败合同禁止启动。

## 资源与归档

资源预检剩余磁盘约252.971GiB，超过20GiB门槛。旧incomplete full upload仍存在，实际占用3.320068359GiB（有重复硬链接，逻辑文件大小合计约6.896531GiB），未删除。结果只包含代码/构建证据、冻结合同、紧凑CSV、网格及失败姿态、日志和核查报告；不含整个临时build树，无周期checkpoint或大量ParaView帧。输出大小见OUTPUT_SIZE.json，小于1GiB。

报告和SHA256结果位于/workspace/hemocell_restore/results/rbc_stage1_20260915_165632，随后下载到/home/lzy/projects/compre_output/rbc_stage1/20260915_165632。远端原始报告/manifest保持不可变；WSL下载SHA256与独立静态几何复核记录见LOCAL_TRANSFER_VERIFICATION.json、LOCAL_GEOMETRY_AUDIT.json与LOCAL_SHA256SUMS。

## 保留的限制与后续

HUMAN_RBC_REVIEW=PENDING。SMALL_PACK_IS_NOT_HCT_VALIDATION=YES。mouse membrane mechanics与interior viscosity为UNVERIFIED；实验mouse RBC标定、RBC_GPU_CORRECTNESS、RBC Stage2、MICROBUBBLE_WITH_RBC均PENDING；READY_FOR_MICROBUBBLE=NO。microbubble、adhesion、ultrasound全部OFF。

NEXT_STEP=USER_REVIEW_BEFORE_RBC_STAGE2；这里首先需要审阅并解决Stage1失败，绝不表示已允许进入Stage2。本任务没有自动运行任何Hct群体、GPU IBM、GPU particle mechanics或微泡任务。
