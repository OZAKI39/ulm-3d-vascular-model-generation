# 当前自适应盒体：2 mm STL

侧壁及底板均为 **2 mm**。内腔尺寸保持为 78.17 × 70.53 × 29.01 mm，盒体外尺寸为 **82.17 × 74.53 × 31.01 mm**。

中央血管和四根延长段直接使用既有选中模型，未重新拟合或选路。内壁位置保持，侧壁向外侧减薄、底板从底面减薄。由于外壁退回，墙外可见接口增加至 I1 7.00 mm、O1 7.00 mm、O2 7.00 mm、O3 7.00 mm，端口总长度没有改变。

- `final_abs_casting_mold_2mm.stl`：包含血管、端口和盒体的一体 STL，保留铸造坐标。
- `final_abs_casting_mold_2mm_oriented.stl`：相同模型，按已选旋转姿态重新居中、底面落在 Z=0，供导入切片器。
- `intact_open_top_box_2mm.stl`：单独盒体。
- `geometry_summary.json`、CSV：实际 STL 厚度、端口穿墙连接、接口和材料保留检查。

模型为单连通、闭合且流形的实体，血管材料损失 0.00000121 mm³。各端口穿墙检查和接口检查均通过。此前 ROI 未达到 60 × 50 mm 容差的限制仍保留。

本次仅重新生成 STL，没有针对新厚度切片。旧目录中的 3MF/G-code 对应原 3 mm 版本，不是本次 2 mm STL 的切片结果。

复现：`./run_python.sh s1-7_support_aware_roi_60x50.py --stage box-stl --box-thickness-mm 2`。
