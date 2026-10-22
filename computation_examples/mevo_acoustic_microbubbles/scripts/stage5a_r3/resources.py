"""Low-frequency userspace resource observations, measured monitoring cost."""

import subprocess
import threading
import time
import psutil


class ResourceMonitor:
    def __init__(self, workers=1):
        self.workers = workers
        self.rows = []
        self.overhead_s = 0.0
        self.stop = threading.Event()
        self.process = psutil.Process()
        self.thread = threading.Thread(target=self.collect, daemon=True)

    def collect(self):
        self.process.cpu_percent()
        psutil.cpu_percent(percpu=True)
        while not self.stop.is_set():
            start = time.perf_counter()
            cores = psutil.cpu_percent(percpu=True)
            row = dict(
                monotonic_s=start,
                rss_bytes=self.process.memory_info().rss,
                process_cpu_percent=self.process.cpu_percent(),
                per_core_percent=cores,
                host_cpu_percent=sum(cores) / max(1, len(cores)),
                system_memory_used_bytes=psutil.virtual_memory().used,
                configured_component_workers=self.workers,
                process_threads=self.process.num_threads(),
                thread_cpu_times=[
                    dict(id=t.id, user_s=t.user_time, system_s=t.system_time) for t in self.process.threads()
                ],
            )
            row["cpu_percent"] = row["process_cpu_percent"]
            try:
                values = (
                    subprocess.check_output(
                        [
                            "nvidia-smi",
                            "--query-gpu=utilization.gpu,memory.used",
                            "--format=csv,noheader,nounits",
                        ],
                        text=True,
                        timeout=2,
                    )
                    .strip()
                    .split(",")
                )
                row.update(gpu_percent=float(values[0]), gpu_memory_MiB=float(values[1]))
            except Exception as exc:
                row["gpu_observation_error"] = type(exc).__name__
            self.overhead_s += time.perf_counter() - start
            self.rows.append(row)
            self.stop.wait(0.5)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join(timeout=3)
