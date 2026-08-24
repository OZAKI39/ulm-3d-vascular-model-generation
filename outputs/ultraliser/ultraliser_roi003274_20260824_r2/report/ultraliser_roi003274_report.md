# Ultraliser ROI 003274 report

1. **编译/执行**：成功；commit `3e4b0eee685adbf513e40720a68fd92e66a34b44`。
2. **正式命令**：`wsl -d Ubuntu -- /mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/Ultraliser/build-wsl/bin/ultraVessMorpho2Mesh --morphology /mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/ultraliser/ultraliser_roi003274_20260824_r2/input/roi_core.h5 --output-directory /mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/ultraliser/ultraliser_roi003274_20260824_r2/final/raw_ultraliser --prefix roi003274_final --scaled-resolution --voxels-per-micron 6 --solid --voxelization-axis xyz --packing-algorithm polylines-with-spheres --isosurface-technique dmc --adaptive-optimization --optimization-iterations 5 --smooth-iterations 5 --laplacian-iterations 10 --export-stl-mesh --stats --threads 8`
3. **Packing algorithm**：`polylines-with-spheres`；paper vascular workflow and source POLYLINE_SPHERE_PACKING branch; help default is polylines。
4. **voxels_per_micron**：`6.0`；由最小直径目标和 150,000,000 voxel 上限自动确定。
5. **正式运行耗时**：`40.0313848000369` s。
6. **Surface topology**：watertight=True，components=1，self_intersections=0，nonmanifold_edges=0。
7. **Radius P95 absolute relative error**：`0.11948067655651681`。
8. **Junction 锯齿人工核验**：The former hard merge corner and visible interface at the selected worst junction are absent; the bifurcation is continuous and the saw-tooth artifact is materially improved. Mild voxel-scale surface corrugation remains.。
9. **相对旧模型**：Ultraliser improves junction continuity, but is not clearly superior overall: its mesh has mild voxel-scale corrugation and all 40 successful branch-local samples overestimate source radius (median +10.15%, P95 absolute relative error 11.95%, maximum +12.28%). This systematic dilation is not acceptable for a formal CFD lumen backend without a validated correction.。
10. **建议**：`ULTRALISER_GEOMETRY_NOT_ACCEPTABLE`。

## Compatibility note

本地 Ultraliser vascular reader 仅接受 H5/VMV，和论文声称的 vascular SWC 直读不一致。`roi_core.swc` 是保持 µm 坐标、半径及 parent-child 的 canonical 输入；正式 executable 通过 `roi_core.h5` 官方 vascular schema 读取同一几何，未修改 Ultraliser C++ 算法。
