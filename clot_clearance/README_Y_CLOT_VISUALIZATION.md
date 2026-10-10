# Y 血栓默认科学可视化

正式默认入口为 `scripts/render_y_clot.py`，配置在 `references/y_clot_default_render.json`。
默认仅使用点阵：原始材料节点为胭脂红，固定锚点为浅青灰，真实游离节点为深灰，存活键为细淡线。
旧 medical shell 程序及结果保留为历史可选输出，不是此入口的默认方式。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python scripts/render_y_clot.py
```

加 `--stills` 只生成静态图。此入口只读取已接受数据，不启动 LAMMPS 或 CFD。
局部近景与完整 Y 管并列；位置、位移均使用原始数值，点径是屏幕符号。
每帧给出绝对物理时间、播放倍率及数据来源清单。
输运阶段附着主体明确标注为 `FROZEN REFERENCE SNAPSHOT AFTER 2.33 MS`。

本轮结果及审计：`results/y_clot_point_cloud_and_early_flow_coupling_v1/`、
`verification/y_clot_point_cloud_and_early_flow_coupling_v1/`。
