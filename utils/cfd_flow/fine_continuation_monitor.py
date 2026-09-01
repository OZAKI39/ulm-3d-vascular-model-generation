"""Low-overhead watchdog for the research-only Fine postfix continuation.

The watchdog observes the already-running Musubi process.  It never launches,
signals, or changes that process.  A complete restart is hashed only once and
cached in the current monitor snapshot.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


RUN_RELATIVE = Path(
    "outputs/cfd_flow/"
    "healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901"
)
RUNTIME_WSL = "/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix"
SEGMENT_WSL = f"{RUNTIME_WSL}/segments/segment_0005000_to_1011895"
EXPECTED_BINARY_WSL = (
    "/home/lzy/apes-worktrees/musubi_adaptive_target_fixed_20260901/build/musubi"
)
EXPECTED_BINARY_SHA256 = "1491e7ded30ed15158bd8ba2812d414ad9fa31cbc0128a26d2b1eba779773caa"
EXPECTED_LAUNCHER_WSL = (
    "/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/"
    f"{RUN_RELATIVE.as_posix()}/run_fine_postfix_long.sh"
)
EXPECTED_RANKS = 4
EXPECTED_PAYLOAD_BYTES = 400_949 * 19 * 8
DT_S = 1.2063526530710723e-9
TARGET_ITERATION = 1_011_895
CHECKPOINTS = (202_379, 404_758, 607_137, 809_516, TARGET_ITERATION)
TERMINAL_STATUSES = {
    "CLEAN_COMPLETION",
    "RECOVERABLE_OPERATIONAL_ERROR",
    "SCIENTIFIC_FAILURE",
    "PROCESS_IDENTITY_MISMATCH",
}

_WSL_PROBE = r"""
import datetime
import json
import os
import pathlib
import re
import shutil
import sys

runtime = pathlib.Path(sys.argv[1])
segment = pathlib.Path(sys.argv[2])
expected_binary = sys.argv[3]
expected_launcher = sys.argv[4]
expected_payload_bytes = int(sys.argv[5])
expected_ranks = int(sys.argv[6])

def cmdline(pid):
    try:
        return pathlib.Path(f"/proc/{pid}/cmdline").read_bytes().rstrip(b"\0").split(b"\0")
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return []

def proc_cwd(pid):
    try:
        return os.path.realpath(f"/proc/{pid}/cwd")
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return None

def proc_stat(pid):
    try:
        raw = pathlib.Path(f"/proc/{pid}/stat").read_text()
        fields = raw[raw.rfind(")") + 2:].split()
        return {"state": fields[0], "ticks": int(fields[11]) + int(fields[12])}
    except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError, IndexError):
        return {"state": None, "ticks": 0}

pid_path = runtime / "long_run_pid.txt"
try:
    launcher_pid = int(pid_path.read_text().strip())
except (FileNotFoundError, ValueError):
    launcher_pid = None
launcher_args = cmdline(launcher_pid) if launcher_pid is not None else []
launcher_cwd = proc_cwd(launcher_pid) if launcher_pid is not None else None
launcher_alive = bool(launcher_args)
launcher_text = " ".join(value.decode(errors="replace") for value in launcher_args)
launcher_match = (
    launcher_alive
    and expected_launcher in launcher_text
    and launcher_cwd == str(segment)
)

ranks = []
for entry in pathlib.Path("/proc").iterdir():
    if not entry.name.isdigit():
        continue
    pid = int(entry.name)
    args = cmdline(pid)
    if not args:
        continue
    executable = args[0].decode(errors="replace")
    if executable != expected_binary:
        continue
    cwd = proc_cwd(pid)
    if cwd != str(segment):
        continue
    stat = proc_stat(pid)
    ranks.append(
        {
            "pid": pid,
            "cwd": cwd,
            "command": " ".join(value.decode(errors="replace") for value in args),
            "state": stat["state"],
            "cpu_ticks": stat["ticks"],
        }
    )
ranks.sort(key=lambda item: item["pid"])

def file_info(path):
    if not path.is_file():
        return {"path": str(path), "exists": False, "size_bytes": None, "mtime": None}
    stat = path.stat()
    return {
        "path": str(path),
        "exists": True,
        "size_bytes": stat.st_size,
        "mtime": datetime.datetime.fromtimestamp(
            stat.st_mtime, datetime.timezone.utc
        ).isoformat(),
    }

restart_candidates = []
for directory in (runtime / "restart", segment / "restart"):
    if not directory.is_dir():
        continue
    for header in directory.glob("*_header_*.lua"):
        match = re.search(r"_header_(\d+)\.lua$", header.name)
        if match is None:
            continue
        iteration = int(match.group(1))
        payloads = [
            path
            for path in directory.glob(f"*_{iteration}.lsb")
            if "header" not in path.name and "lastHeader" not in path.name
        ]
        if len(payloads) != 1:
            continue
        payload = payloads[0]
        try:
            header_text = header.read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            continue
        header_iteration_match = re.search(r"iter\s*=\s*(\d+)", header_text)
        header_iteration = (
            int(header_iteration_match.group(1)) if header_iteration_match else None
        )
        complete = (
            header.stat().st_size > 100
            and payload.stat().st_size == expected_payload_bytes
            and header_iteration == iteration
        )
        if complete:
            restart_candidates.append(
                {
                    "iteration": iteration,
                    "header_path": str(header),
                    "path": str(payload),
                    "size_bytes": payload.stat().st_size,
                    "mtime": datetime.datetime.fromtimestamp(
                        payload.stat().st_mtime, datetime.timezone.utc
                    ).isoformat(),
                }
            )
restart_candidates.sort(key=lambda item: item["iteration"])
latest_restart = restart_candidates[-1] if restart_candidates else None

stdout_path = segment / "musubi_stdout.log"
stderr_path = segment / "musubi_stderr.log"
stdout_info = file_info(stdout_path)
stderr_info = file_info(stderr_path)
stdout_text = stdout_path.read_text(encoding="utf-8", errors="replace") if stdout_path.is_file() else ""
stderr_tail = ""
if stderr_path.is_file():
    with stderr_path.open("rb") as stream:
        stream.seek(max(0, stderr_path.stat().st_size - 524288))
        stderr_tail = stream.read().decode("utf-8", errors="replace")
main_text = stdout_text.split("Starting Musubi MAIN loop", 1)[-1] if "Starting Musubi MAIN loop" in stdout_text else ""
adaptive_iterations = [int(value) for value in re.findall(r"ADAPTIVE_FLUX_PRESSURE iter=(\d+)", main_text)]
reported_iterations = [int(value) for value in re.findall(r"(?m)^\s*iterations:\s*(\d+)\s*$", main_text)]
observed = []
if adaptive_iterations:
    observed.append(adaptive_iterations[-1])
if reported_iterations:
    observed.append(reported_iterations[-1])
if latest_restart is not None:
    observed.append(latest_restart["iteration"])
current_iteration = max(observed) if observed else None
sim_times = re.findall(r"(?m)^\s*simTime\s*:\s*([+\-0-9.Ee]+)\s*$", main_text)
physical_time = float(sim_times[-1]) if sim_times else None

fatal_patterns = {
    "nan": r"(?i)(?<![A-Za-z])nan(?![A-Za-z])",
    "segmentation_fault": r"(?i)segmentation fault",
    "floating_point_exception": r"(?i)floating point exception",
    "mpi_abort": r"(?i)MPI_ABORT",
    "out_of_memory": r"(?i)out of memory|oom-kill",
    "negative_pdf": r"(?i)negative pdf|non-positive pdf",
    "solver_abort": r"(?i)(?:solver|musubi).*\babort(?:ed|ing)?\b",
}
fatal_matches = []
fatal_text = main_text + "\n" + stderr_tail
for label, pattern in fatal_patterns.items():
    if re.search(pattern, fatal_text):
        fatal_matches.append(label)

mem_available = None
for line in pathlib.Path("/proc/meminfo").read_text().splitlines():
    if line.startswith("MemAvailable:"):
        mem_available = int(line.split()[1]) * 1024
        break
pressure = {}
for kind in ("cpu", "memory", "io"):
    path = pathlib.Path(f"/proc/pressure/{kind}")
    if not path.is_file():
        pressure[kind] = None
        continue
    first = path.read_text().splitlines()[0]
    match = re.search(r"avg10=([0-9.]+)", first)
    pressure[kind] = float(match.group(1)) if match else None

runtime_status_path = runtime / "long_run_status.txt"
runtime_status = runtime_status_path.read_text(errors="replace").strip() if runtime_status_path.is_file() else None
successful_run = "SUCCESSFUL run" in stdout_text
payload = {
    "launcher_pid": launcher_pid,
    "launcher_alive": launcher_alive,
    "launcher_command": launcher_text,
    "launcher_cwd": launcher_cwd,
    "launcher_match": launcher_match,
    "ranks": ranks,
    "expected_rank_count": expected_ranks,
    "command_match": launcher_match and len(ranks) == expected_ranks,
    "cpu_total_ticks": sum(item["cpu_ticks"] for item in ranks),
    "clock_ticks_per_second": os.sysconf(os.sysconf_names["SC_CLK_TCK"]),
    "system_uptime_s": float(pathlib.Path("/proc/uptime").read_text().split()[0]),
    "load_average": list(os.getloadavg()),
    "available_memory_bytes": mem_available,
    "free_disk_bytes": shutil.disk_usage(runtime).free,
    "pressure_avg10": pressure,
    "latest_restart": latest_restart,
    "complete_restart_count": len(restart_candidates),
    "stdout": stdout_info,
    "stderr": stderr_info,
    "current_iteration": current_iteration,
    "physical_time_s": physical_time,
    "fatal_matches": fatal_matches,
    "runtime_status": runtime_status,
    "successful_run": successful_run,
}
print(json.dumps(payload, sort_keys=True))
"""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def _run_probe() -> dict[str, Any]:
    process = subprocess.run(
        [
            "wsl.exe",
            "-d",
            "Ubuntu",
            "--",
            "python3",
            "-c",
            _WSL_PROBE,
            RUNTIME_WSL,
            SEGMENT_WSL,
            EXPECTED_BINARY_WSL,
            EXPECTED_LAUNCHER_WSL,
            str(EXPECTED_PAYLOAD_BYTES),
            str(EXPECTED_RANKS),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    if process.returncode != 0:
        raise RuntimeError(process.stderr.strip() or f"WSL probe returned {process.returncode}")
    return json.loads(process.stdout)


def _sha256_wsl(path: str) -> str:
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "sha256sum", path],
        check=False,
        capture_output=True,
        text=True,
        timeout=180,
    )
    if process.returncode != 0:
        raise RuntimeError(process.stderr.strip() or "restart sha256 failed")
    digest = process.stdout.split()[0]
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise RuntimeError("restart sha256 output was malformed")
    return digest


def _cpu_activity(probe: dict[str, Any], previous: dict[str, Any]) -> dict[str, Any]:
    previous_probe = previous.get("_probe", {})
    ticks = int(probe.get("cpu_total_ticks", 0))
    previous_ticks = previous_probe.get("cpu_total_ticks")
    uptime = float(probe.get("system_uptime_s", 0.0))
    previous_uptime = previous_probe.get("system_uptime_s")
    aggregate_percent = None
    if previous_ticks is not None and previous_uptime is not None and uptime > float(previous_uptime):
        cpu_seconds = (ticks - int(previous_ticks)) / float(probe["clock_ticks_per_second"])
        aggregate_percent = 100.0 * cpu_seconds / (uptime - float(previous_uptime))
    states = [item.get("state") for item in probe.get("ranks", [])]
    active = bool(states) and all(state not in {None, "Z", "X"} for state in states)
    return {
        "aggregate_percent": aggregate_percent,
        "rank_count": len(states),
        "rank_states": states,
        "active": active,
    }


def _classify(probe: dict[str, Any]) -> str:
    fatal = bool(probe.get("fatal_matches"))
    ranks_alive = len(probe.get("ranks", [])) == EXPECTED_RANKS
    if fatal:
        return "SCIENTIFIC_FAILURE"
    if ranks_alive and probe.get("command_match"):
        return "HEALTHY_RUNNING"
    if probe.get("launcher_alive") or ranks_alive:
        return "PROCESS_IDENTITY_MISMATCH"
    latest = probe.get("latest_restart") or {}
    complete_target = int(latest.get("iteration", -1)) == TARGET_ITERATION
    if probe.get("runtime_status") == "PASS" and complete_target and probe.get("successful_run"):
        return "CLEAN_COMPLETION"
    return "RECOVERABLE_OPERATIONAL_ERROR"


def _next_checkpoint(latest_iteration: int | None) -> int | None:
    current = -1 if latest_iteration is None else latest_iteration
    return next((value for value in CHECKPOINTS if value > current), None)


def restart_pair_complete(
    *, header_text: str, header_size: int, payload_size: int, iteration: int
) -> bool:
    """Mirror the watchdog's deliberately cheap restart completeness gate."""

    match = re.search(r"iter\s*=\s*(\d+)", header_text)
    return (
        header_size > 100
        and payload_size == EXPECTED_PAYLOAD_BYTES
        and match is not None
        and int(match.group(1)) == iteration
    )


def build_record(probe: dict[str, Any], previous: dict[str, Any]) -> dict[str, Any]:
    """Build one monitor record and retain the restart hash cache."""

    hashes = dict(previous.get("_restart_sha256_by_iteration", {}))
    latest = probe.get("latest_restart") or {}
    latest_iteration = latest.get("iteration")
    if latest_iteration is not None:
        key = str(int(latest_iteration))
        if key not in hashes:
            hashes[key] = _sha256_wsl(str(latest["path"]))
        latest_hash = hashes[key]
    else:
        latest_hash = None
    current_iteration = probe.get("current_iteration")
    physical_time = probe.get("physical_time_s")
    if physical_time is None and current_iteration is not None:
        physical_time = int(current_iteration) * DT_S
    previous_probe = previous.get("_probe", {})
    record = {
        "timestamp": _utc_now(),
        "pid": probe.get("launcher_pid"),
        "process_alive": bool(probe.get("launcher_alive"))
        and len(probe.get("ranks", [])) == EXPECTED_RANKS,
        "command_match": bool(probe.get("command_match")),
        "current_iteration": current_iteration,
        "current_iteration_observation": "LATEST_EMITTED_OR_COMPLETE_CHECKPOINT",
        "physical_time_s": physical_time,
        "latest_restart_iteration": latest_iteration,
        "latest_restart_path": latest.get("path"),
        "latest_restart_header_path": latest.get("header_path"),
        "latest_restart_sha256": latest_hash,
        "latest_restart_size_bytes": latest.get("size_bytes"),
        "stdout_last_update": (probe.get("stdout") or {}).get("mtime"),
        "stderr_last_update": (probe.get("stderr") or {}).get("mtime"),
        "stdout_updated_since_previous": (probe.get("stdout") or {}).get("mtime")
        != (previous_probe.get("stdout") or {}).get("mtime"),
        "stderr_updated_since_previous": (probe.get("stderr") or {}).get("mtime")
        != (previous_probe.get("stderr") or {}).get("mtime"),
        "cpu_activity": _cpu_activity(probe, previous),
        "available_memory_bytes": probe.get("available_memory_bytes"),
        "free_disk_bytes": probe.get("free_disk_bytes"),
        "load_average": probe.get("load_average"),
        "resource_pressure_avg10": probe.get("pressure_avg10"),
        "fatal_pattern_found": bool(probe.get("fatal_matches")),
        "fatal_patterns": probe.get("fatal_matches", []),
        "next_expected_checkpoint": _next_checkpoint(latest_iteration),
        "status": _classify(probe),
        "runtime_status": probe.get("runtime_status"),
        "rank_pids": [item["pid"] for item in probe.get("ranks", [])],
        "expected_binary": EXPECTED_BINARY_WSL,
        "expected_binary_sha256": EXPECTED_BINARY_SHA256,
        "_restart_sha256_by_iteration": hashes,
        "_probe": probe,
    }
    return record


def _public_record(record: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def _git_head(project_root: Path) -> str | None:
    process = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=project_root,
        check=False,
        capture_output=True,
        text=True,
        timeout=10,
    )
    return process.stdout.strip() if process.returncode == 0 else None


def _write_oracle_state(project_root: Path, record: dict[str, Any]) -> None:
    qc = project_root / RUN_RELATIVE / "qc"
    oracle_path = qc / "adaptive_active_population_cb_equivalence.json"
    oracle = _read_json(oracle_path)
    if not oracle:
        oracle = {
            "status": "CB_ORACLE_DEFERRED_UNTIL_FINE_COMPLETION",
            "reason": "Fine is the primary task; the oracle has not been launched.",
            "solver_calls": 0,
        }
    fine_status = str(record["status"])
    if fine_status == "HEALTHY_RUNNING":
        next_action = "CONTINUE_120_SECOND_FINE_MONITORING"
    elif fine_status == "CLEAN_COMPLETION":
        next_action = "FREEZE_FINAL_RESTART_AND_RUN_FINE_STEADY_AUDITOR"
    elif fine_status == "RECOVERABLE_OPERATIONAL_ERROR":
        next_action = "VERIFY_COMPATIBILITY_BEFORE_ANY_CHECKPOINT_RESUME"
    else:
        next_action = "STOP_AND_PRESERVE_EVIDENCE"
    _write_json_atomic(
        qc / "monitor_and_cb_oracle_state.json",
        {
            "timestamp": record["timestamp"],
            "git_head": _git_head(project_root),
            "fine": _public_record(record),
            "cb_oracle": oracle,
            "next_action": next_action,
        },
    )


def poll_once(project_root: Path) -> dict[str, Any]:
    qc = project_root / RUN_RELATIVE / "fine_postfix" / "qc"
    current_path = qc / "fine_continuation_monitor.json"
    history_path = qc / "fine_continuation_monitor.jsonl"
    previous = _read_json(current_path)
    probe = _run_probe()
    record = build_record(probe, previous)
    _write_json_atomic(current_path, record)
    with history_path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(_public_record(record), sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    _write_oracle_state(project_root, record)
    return record


def monitor(project_root: Path, interval_s: float, once: bool) -> int:
    while True:
        try:
            record = poll_once(project_root)
        except Exception as exc:  # keep a transient WSL probe failure observable
            qc = project_root / RUN_RELATIVE / "fine_postfix" / "qc"
            failure = {
                "timestamp": _utc_now(),
                "status": "MONITOR_PROBE_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            }
            _write_json_atomic(qc / "fine_continuation_monitor.json", failure)
            with (qc / "fine_continuation_monitor.jsonl").open(
                "a", encoding="utf-8", newline="\n"
            ) as stream:
                stream.write(json.dumps(failure, sort_keys=True) + "\n")
            if once:
                raise
        else:
            print(json.dumps(_public_record(record), sort_keys=True), flush=True)
            if once or record["status"] in TERMINAL_STATUSES:
                return 0
        time.sleep(interval_s)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--interval", type=float, default=120.0)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    if args.interval <= 0:
        parser.error("--interval must be positive")
    return monitor(args.project_root.resolve(), args.interval, args.once)


if __name__ == "__main__":
    raise SystemExit(main())
