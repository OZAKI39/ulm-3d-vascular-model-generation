# F10S：Y近场零点与PDMS计算适用性

所有Y结论为 **HOMOGENEOUS_WATER_PROXY_SPECIFIC**，几何偏移仅为探索性分析，不能当作用户装配误差。

固定1500 Pa-equivalent、2.5 MHz，Y−源距中心4.5 mm，理论轴上第二零点在4.504166667 mm。源距偏移−25/−5/0/+5/+25 µm，团心压力分别96.808/30.335/13.778/2.753/68.614 Pa；基线Y/X幅值比0.00494616（约1/202.18）。同样的±5 µm合法团心位移也跨越明显压力变化；一次Fy在零点两侧反号。单源公共相位0/90/180°不改变压力模和周期平均力，已检查。

2.475/2.5/2.525 MHz仅场级变化得到174.756/13.778/202.258 Pa。非2.5 MHz没有调用旧APECSS力或运动，保存NOT_RUN_FREQUENCY_OUTSIDE_APECSS。因此Y低效高度依赖波长和此几何的均匀水近场零点，不能外推真实PDMS。

参考材料：Guangyao Xu et al., *Acoustic Characterization of Polydimethylsiloxane for Microscale Acoustofluidics*, Physical Review Applied 13, 054069 (2020), DOI [10.1103/PhysRevApplied.13.054069](https://doi.org/10.1103/PhysRevApplied.13.054069)。Tables II/III的Sylgard184 10:1、85°C固化、25°C测量参考组，测量频带3–7 MHz；原F10M保存论文来源、摘要参数、PDF SHA和表格目视审查。F10S另核查作者/DTU公开记录。论文含Kramers–Kronig拟合，允许作为新假设研究外推，但本轮没有把2.5 MHz外推材料静默准入。

预检保留32.5×9×9 mm块、1.6 mm通道、两个1.5 mm顶部水孔、四7 mm圆片、X封闭几何。域含条件性1 mm外部水包层，BLI额外36网格支持余量，域外20网格PML；参考低频声速仅作条件性资源规划下界，不作为已验证2.5 MHz材料输入。水波长0.6 mm；CFL≤0.2；PDMS等效流体不包含完整剪切波，真实端封材料/厚度与外部耦合仍未校准。

| PPW | 网格间距/µm | 含PML格数 | 估计峰值/GiB | 单源估计/s |
|---|---:|---|---:|---:|
| 6 | 68.112 | 625×275×275 | 60.022 | 1955.8 |
| 8 | 51.084 | 792×330×330 | 109.526 | 4481.9 |

估算按F10L实测三源373.706 s、15100 MiB、378×192×192、3504步，按体素/时间步放缩，显存×1.2、耗时×1.5异质余量。是估计而非本轮实测；此源支持与均匀网格方案并非所有可能方法的最小成本。两候选都超20 GiB安全预算、24 GiB硬件及本轮120 s诊断预算。资源和材料准入未满足，**PDMS_Y_FIELD_NOT_ADMITTED / PDMS_FIELD_NOT_RUN**；没有求解多材料压力场，不伪造对照图。

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED；training_ready=false。

