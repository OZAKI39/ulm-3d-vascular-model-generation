# 刚体运动与真正形变的分离

使用固定ID Kabsch拟合 proper rotation，不允许反射；依次报告原始位移、去平移、去平移/旋转残差，再拟合一个统一尺度，得到 similarity RMS。后者不把各向同性收缩误认为各向异性形变。全局F采用极分解F=RU；U的无迹部分、非对角部分和非仿射残差互补，不能用无效F概括严重重排。

单元测试含纯平移、刚体旋转、统一缩放、已知剪切、非仿射运动、秩退化、近简并特征值和投影假接触。算法参照 [SciPy Kabsch](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html) 与 [极分解](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.polar.html)。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
