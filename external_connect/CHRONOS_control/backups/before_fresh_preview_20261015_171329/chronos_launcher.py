"""Chinese desktop launcher for the local Chronos 1.4 tools."""
import ipaddress
import json
import os
from pathlib import Path
import subprocess
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk
import webbrowser

from chronos_common import BASE, DEFAULT_CONFIG, load_config, video_options, write_json


def main():
    root = tk.Tk()
    root.title("Chronos 1.4 · 本地采集")
    root.geometry("700x575")
    root.minsize(680, 560)
    box = ttk.Frame(root, padding=22)
    box.pack(fill="both", expand=True)
    ttk.Label(box, text="Chronos 1.4", font=("Microsoft YaHei UI", 22, "bold")).pack(anchor="w")
    ttk.Label(box, text="实时预览 · 截图 · 接收帧率 · 标定与测量", font=("Microsoft YaHei UI", 10)).pack(anchor="w", pady=(0, 18))
    settings = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8-sig"))
    ip = tk.StringVar(value=settings["camera_ip"])
    transport = tk.StringVar(value=settings["rtsp_transport"])
    mode_labels = {'low':'低延迟（推荐）', 'compatibility':'兼容模式（原缓存）'}
    mode = tk.StringVar(value=mode_labels[settings.get('latency_mode', 'low')])
    row = ttk.Frame(box)
    row.pack(fill="x")
    ttk.Label(row, text="相机 IP").pack(side="left")
    ttk.Entry(row, textvariable=ip, width=20).pack(side="left", padx=10)
    ttk.Label(row, text="传输").pack(side="left")
    ttk.Combobox(row, textvariable=transport, values=["udp", "tcp"], state="readonly", width=6).pack(side="left", padx=10)
    profile_row = ttk.Frame(box)
    profile_row.pack(fill='x', pady=(10, 0))
    ttk.Label(profile_row, text='预览模式').pack(side='left')
    ttk.Combobox(profile_row, textvariable=mode, values=list(mode_labels.values()),
                 state='readonly', width=23).pack(side='left', padx=10)
    jobs = []
    status = tk.StringVar(value="就绪。当前相机已验证使用 UDP；每次只打开一个视频接收程序。")

    def save_settings():
        try:
            settings["camera_ip"] = str(ipaddress.IPv4Address(ip.get().strip()))
            settings["rtsp_transport"] = transport.get()
            settings['latency_mode'] = next(key for key,value in mode_labels.items() if value == mode.get())
            write_json(DEFAULT_CONFIG, settings)
            return load_config()
        except Exception as exc:
            messagebox.showerror("设置有误", str(exc), parent=root)
            return None

    def run(script, arguments=(), video=False):
        config = save_settings()
        if config is None:
            return
        if video and any(is_video and process.poll() is None for process, handle, is_video, name in jobs):
            messagebox.showinfo("视频正在运行", "请先在预览窗口按 Q 退出，或等待稳定性测试结束。", parent=root)
            return
        logs = Path(config["output_dir"]) / "launcher_logs"
        logs.mkdir(parents=True, exist_ok=True)
        from datetime import datetime
        log = logs / f"{Path(script).stem}_{datetime.now():%Y%m%d_%H%M%S_%f}.log"
        handle = log.open("w", encoding="utf-8")
        environment = os.environ.copy()
        environment["OPENCV_FFMPEG_CAPTURE_OPTIONS"] = video_options(config)
        environment["PYTHONIOENCODING"] = "utf-8"
        try:
            process = subprocess.Popen([sys.executable, "-u", str(BASE / script), *map(str, arguments)],
                                       cwd=BASE, env=environment, stdout=handle, stderr=subprocess.STDOUT,
                                       creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        except Exception as exc:
            handle.close()
            messagebox.showerror("启动失败", str(exc), parent=root)
            return
        jobs.append((process, handle, video, script))
        status.set(f"已启动 {script}。运行日志保存在 data\\launcher_logs。")

    def calibrate():
        image = filedialog.askopenfilename(title="选择测微尺原始截图", initialdir=BASE / "data", filetypes=[("图像", "*.png *.jpg *.tif")])
        if not image:
            return
        distance = simpledialog.askfloat("已知长度", "准备点击的两个刻度之间相距多少微米？", minvalue=0.000001, parent=root)
        if distance is not None:
            note = simpledialog.askstring("记录光学配置", "填写物镜倍率及其他需要记录的设置：", parent=root) or ""
            run("calibrate_scale.py", ["--image", image, "--distance-um", distance, "--note", note])

    def analyze():
        directory = filedialog.askdirectory(title="选择连续截图的 images 文件夹", initialdir=BASE / "data")
        if not directory:
            return
        roi = simpledialog.askstring("测量区域", "填写区域 x y 宽 高（原图像素）；留空使用整图。\n程序测量区域内最大的完整深色轮廓，请检查结果叠加图。", parent=root)
        if roi is None:
            return
        arguments = ["--images", directory]
        if roi.strip():
            try:
                values = [int(v) for v in roi.replace(",", " ").split()]
                if len(values) != 4:
                    raise ValueError()
            except ValueError:
                messagebox.showerror("区域有误", "请输入四个整数，例如 500 300 200 200。", parent=root)
                return
            arguments += ["--roi", *values]
        if (BASE / "calibration.json").exists() and messagebox.askyesno("使用标定", "当前物镜、分辨率、ROI 和光路是否与 calibration.json 标定时完全一致？", parent=root):
            arguments += ["--calibration", BASE / "calibration.json"]
        run("analyze_stability.py", arguments)

    buttons = ttk.Frame(box)
    buttons.pack(fill="x", pady=18)
    actions = [
        ("打开实时画面", lambda: run("test_chronos.py", video=True)),
        ("连接诊断", lambda: run("check_connection.py")),
        ("5 分钟稳定性测试", lambda: run("test_chronos.py", ["--headless", "--seconds", 300, "--save-every", 60, "--max-images", 5], video=True)),
        ("采样 60 秒（每秒一图）", lambda: run("test_chronos.py", ["--seconds", 60, "--save-every", 1, "--max-images", 60], video=True)),
        ("打开相机网页", lambda: save_settings() and webbrowser.open(f"http://{ip.get().strip()}")),
        ("打开数据文件夹", lambda: open_data()),
        ("测微尺标定", calibrate),
        ("分析截图中的轮廓", analyze),
        ("显示延迟测量时钟", lambda: run("measure_latency.py", ["--clock"])),
        ("打开使用说明", lambda: os.startfile(BASE / "使用说明.txt")),
        ("使用 VLC 预览", lambda: run("open_vlc.py", video=True)),
        ("保存相机设置", save_settings),
    ]
    for i, (label, action) in enumerate(actions):
        ttk.Button(buttons, text=label, command=action).grid(row=i//2, column=i%2, sticky="ew", padx=4, pady=5)
    buttons.columnconfigure(0, weight=1)
    buttons.columnconfigure(1, weight=1)
    ttk.Label(box, text="预览窗口：S 后台保存原图，R 重新连接刷新，Q 退出。切换模式后重新打开预览。", wraplength=650).pack(anchor="w")
    ttk.Label(box, textvariable=status, wraplength=610, foreground="#145b72").pack(anchor="w", pady=10)

    def open_data():
        config = save_settings()
        if config:
            Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)
            os.startfile(config["output_dir"])

    def poll():
        for job in jobs[:]:
            process, handle, is_video, name = job
            if process.poll() is not None:
                handle.close()
                jobs.remove(job)
                status.set(f"{name} 已结束，返回值 {process.returncode}。结果和日志在 data 文件夹。")
        root.after(500, poll)

    def close():
        if any(p.poll() is None for p, _, _, _ in jobs):
            messagebox.showinfo("任务仍在运行", "请先在预览窗口按 Q 退出；等待后台测试结束后再关闭面板。", parent=root)
            return
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", close)
    poll()
    root.mainloop()


if __name__ == "__main__":
    main()
