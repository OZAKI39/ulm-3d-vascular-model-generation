# 复现命令

以下两条命令只读取既有 WSS、端口和管网文件，输出到当前 meeting_question_maps/；不调用求解器或任务管理器。

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 \
  /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B \
  /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/wss_validation_v2/meeting_question_maps/scripts/make_meeting_maps.py \
  > /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/wss_validation_v2/meeting_question_maps/logs/render.log 2>&1

PYTHONDONTWRITEBYTECODE=1 /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B \
  /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/wss_validation_v2/meeting_question_maps/scripts/write_explanation.py
```

中文字体使用 `/mnt/c/Windows/Fonts/msyh.ttc` 和 `msyhbd.ttc`；换环境时可改成已安装的中文字体路径。Matplotlib 缓存固定放在本目录 `.mplconfig/`。PyVista 使用离屏渲染。
