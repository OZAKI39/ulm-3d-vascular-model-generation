"""Read-only Windows USB + VISA diagnosis; saves connection_report.json."""
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from afg_control import BASE, discover, load_config


def main():
    report = {"checked_at_utc": datetime.now(timezone.utc).isoformat(), "python": sys.executable,
              "platform": platform.platform(), "note": "Read-only: no output or parameter changes."}
    if sys.platform == "win32":
        script = r"""
        $ErrorActionPreference = 'Stop'
        Get-PnpDevice -PresentOnly | Where-Object { $_.InstanceId -match 'VID_0699&PID_0353' } |
          ForEach-Object {
            $problem = Get-PnpDeviceProperty -InstanceId $_.InstanceId -KeyName DEVPKEY_Device_ProblemCode
            [PSCustomObject]@{Name=$_.FriendlyName;Status=$_.Status;InstanceId=$_.InstanceId;ProblemCode=$problem.Data}
          } | ConvertTo-Json -Depth 3 -Compress
        """
        try:
            result = subprocess.run(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", script],
                                    capture_output=True, timeout=40, creationflags=subprocess.CREATE_NO_WINDOW)
            text = result.stdout.decode("utf-8", errors="replace").strip()
            report["windows_usb"] = json.loads(text) if result.returncode == 0 and text else {"error": result.stderr.decode(errors="replace")}
        except Exception as e:
            report["windows_usb"] = {"error": str(e)}
    try:
        report["visa_instruments"] = discover(load_config())
    except Exception as e:
        report["visa_error"] = str(e)
    target = BASE / "connection_report.json"
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print("Saved:", target)
    return 1 if "visa_error" in report else 0


if __name__ == "__main__":
    raise SystemExit(main())
