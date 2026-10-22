# RFUS 永久合成试件

`manifest.json` 列出 15 个几何试件；同名 `*_vascular.stl` 为血管，`*_obstacles.stl` 为碰撞或受困障碍。它们由 `vascular_processing/support_rfus_fixtures.py` 生成，测试代码在 `tests/test_s1_8_vascular_support_planner.py` 中。所有试件是合成几何，不是 BG001 生产血管。

真实 Bambu 长水平管试验另外保存在：

`outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner/robust_full_underbody/validation/actual_straight_tube/`

该试验使用 30 mm 水平圆管、双线帘、双线底颈、2 mm 试验底板。远离圆管投影的独立竖柱保证空气间隙层不是整个试件的空层；它不计入下托覆盖。覆盖根据归档 G-code 的实际线宽测量。除从测得足迹中去掉中间 20 mm 的单元负对照外，还生成并实际切片一件中间缺少 20 mm 帘体的稀疏对照，要求根据其归档 G-code 测出缺口。试验不证明帘强度或实物可拆除性。

统一验证入口：

```bash
./run_python.sh tools/verify_support_rfus.py
```

几何测试与真实切片证据检查分开记录；只有前者通过不能替代后者。
