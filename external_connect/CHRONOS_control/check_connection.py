"""Read-only Chronos HTTP/API/RTSP checks; never changes the camera."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import socket
import sys
import urllib.error
import urllib.request

from chronos_common import DEFAULT_CONFIG, load_config, session_dir, utc_now, write_json


def port_check(ip, port, timeout=2):
    try:
        with socket.create_connection((ip, port), timeout):
            return {"reachable": True}
    except OSError as exc:
        return {"reachable": False, "error": str(exc)}


def get_http(ip, path, timeout=3):
    # Camera traffic stays local and bypasses corporate HTTP proxy settings.
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    try:
        with opener.open(f"http://{ip}{path}", timeout=timeout) as response:
            raw = response.read(2_000_000).decode("utf-8", errors="replace")
            try:
                body = json.loads(raw)
            except ValueError:
                body = raw[:500]
            return {"ok": True, "http_status": response.status, "body": body}
    except (OSError, urllib.error.URLError) as exc:
        return {"ok": False, "error": str(exc)}


def diagnose(config):
    ip = config["camera_ip"]
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda port: port_check(ip, port), [80, 554]))
    report = {"time_utc": utc_now(), "camera_ip": ip,
              "http_port_80": results[0], "rtsp_port_554": results[1],
              "web_url": f"http://{ip}", "rtsp_url": f"rtsp://{ip}",
              "note": "An open RTSP port does not prove successful video decoding."}
    if results[0]["reachable"]:
        paths = ["/", "/apidoc/", "/control/p/cameraModel", "/control/p/cameraSerial",
                 "/control/p/cameraApiVersion", "/control/p/state", "/control/p/videoState",
                 "/control/p/resolution", "/control/p/framePeriod", "/control/p/exposurePeriod"]
        with ThreadPoolExecutor(max_workers=3) as pool:
            report["http"] = dict(zip(paths, pool.map(lambda p: get_http(ip, p), paths)))
    try:
        import cv2
        report["opencv_version"] = cv2.__version__
        report["ffmpeg_enabled"] = any("FFMPEG:" in line and "YES" in line for line in cv2.getBuildInformation().splitlines())
    except ImportError:
        report["opencv_version"] = None
        report["ffmpeg_enabled"] = False
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--ip", help="Override camera IP for this run")
    args = parser.parse_args()
    config = load_config(args.config, args.ip)
    report = diagnose(config)
    path = session_dir(config, "connection") / "report.json"
    write_json(path, report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print(f"Report: {path}")
    return 0 if report["http_port_80"]["reachable"] and report["rtsp_port_554"]["reachable"] else 2


if __name__ == "__main__":
    sys.exit(main())
