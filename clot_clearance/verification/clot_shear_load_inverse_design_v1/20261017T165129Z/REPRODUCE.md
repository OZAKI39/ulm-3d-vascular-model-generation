# 当前表面随动剪切复现

已接受目录 `verification/clot_shear_load_inverse_design_v1/20261017T165129Z` 只读。新计算先复制分析目录并设置新的 common.py RUN、Supervisor名称，重新冻结，不覆盖历史结果。

环境：`source /workspace/y_clot/env.sh`；`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`。仅1 rank。沿用 `verification/clot_multifragment_breakup_v1/20261017T110849Z/build/liblammps.so`，SHA 465f8144c6da5c9de3aec4137b9c325f53c26737b199c536fb4a20845511ef77，LD_LIBRARY_PATH加入该build。无需重编译或安装新runtime依赖。构建配方在前序BUILD_AUDIT.json与analysis/clot_multifragment_breakup_v1/build.py。

执行顺序（服务器工程根目录，脚本在本轮analysis目录）：

1. `python prepare.py`：冻结原试样、材料、参考证据、两种场参数、原始面片和独立碎块标准。只复用当前调用链。来源见BASELINE_MANIFEST.json。
2. `python checks.py`：新路径几何/切向/方向/力矩/功率/原始ID、解析梯度和两个微小native步、同态restart force验证。不是新的主断裂工况。
3. `python freeze_sources.py`：主计算前源码、输入、native版本和git SHA冻结。不得随后修改冻结的11个运行/准备文件。
4. `python manage.py uniform`、`python manage.py counter_rotating`：分别Supervisor前台运行，可并行，每例1线程。dt=.25µs、事件最小.03125µs，观察3.6ms。原始结果若存在则拒绝覆盖。没有固定4GiB或1800s限制。遇数值错误保存失败状态和最后接受态，不将失败态混入正式结果。
5. `python analyze.py`：真实存活图、锚可达、尺寸/寿命/内部连接、自然失效网络、只读桥接/谱候选、表面覆盖和能量分析。
6. `python validate_outputs.py`：每例一次逐接受步原始ID损伤与完整连通性重放、FP64载荷回读、材料守恒和账本重建。可用 `--case uniform` 或 `--case counter_rotating` 分别执行，按原始SHA合并证据，不重新动力学。
7. `python surface_coverage_audit.py` 区分非表面内点与整个无受载面片分量，随后 `python figures.py`、`python report.py`：PNG/SVG、20问中文报告及HTML。`--inputs-only`仅先绘制冻结场。所有断裂图来自原生接受态，实际坐标×1。
8. 同步本轮A/O/F到WSL，本地 `python3 update_records.py` 更新科研记录和当前入口，逐字保留历史，写前备份/哈希、diff、读回；上传两个更新证据及入口。人工逐图查看后记录VISUAL_REVIEW.json。
9. `python seal.py`：复核旧依赖、当前源、两服务退出，封存manifest。同步回WSL，`python3 local_readback.py`逐文件SHA及HTML链接回读，回传证明。

固定物理近似：τ=(µ/h)Pₜ(u_outer−v)，µ=.001Pa·s/h25µm。1000Pa是静止面驱动cap，对应25m/s，非实际两例相等剪切/等能量。初始182顶面原节点三角形只有三原边live且非退化才保留；失去面片不删除材料键。每力回调重算当前面积/法向/速度，正1/3权重；速度相关末半步固定点解与native VV一致。2.8ms驱动停止，余下.8ms保留被动拖曳。

新裂面、其他原外表面及无保留面片碎片的水动力未解析。平面涡核心无散不等于变形曲面/截断域满足三维Navier–Stokes。机械标签和适用域诊断必须保留。旧V3有符号映射未修复，没有完整V3/CFD/泡流求解。

输出包括每步标量history、lossless原始ID损伤变化及组件成员、关键FP64坐标/速度/力/κ/d/能量，外功drive/drag分开、当前x×F角冲量、接触激活J及独立R。并行两例的累计算例墙钟不能当成项目经过时间，各自RSS峰值之和仅是并发上界。

原生库及上游构建树保留服务器，交付新任务源码、参数、日志和必要真实轨迹，不迁移无关大文件。local_readback不冒称下载全部依赖。新优化必须依真实20–200碎片或明确近分离证据准入，不以断键数代替。
