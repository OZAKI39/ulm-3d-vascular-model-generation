"""Open the bundled portable VLC with the verified Chronos UDP parameters."""
import subprocess
import argparse
from chronos_common import BASE, load_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=0, help='Stop after N seconds; 0 = manual close')
    args = parser.parse_args()
    if args.seconds < 0:
        parser.error('seconds must be nonnegative')
    config = load_config()
    executable = BASE / "tools" / "vlc-3.0.24" / "vlc.exe"
    if not executable.exists():
        raise FileNotFoundError(f"Portable VLC was not found: {executable}")
    arguments = [str(executable), "--ignore-config", "--no-one-instance", "--no-qt-privacy-ask",
                 "--no-metadata-network-access", "--no-qt-updates-notif", "--network-caching=500",
                 "--no-video-title-show", "--avcodec-hw=none",
                 "--rtsp-tcp" if config["rtsp_transport"] == "tcp" else "--no-rtsp-tcp"]
    if config["rtsp_transport"] == "udp" and config.get("udp_port"):
        arguments.append(f"--rtp-client-port={config['udp_port']}")
    if args.seconds:
        arguments.extend([f'--run-time={args.seconds}', '--play-and-exit'])
    arguments.append(f"rtsp://{config['camera_ip']}")
    return subprocess.call(arguments, cwd=executable.parent)


if __name__ == "__main__":
    raise SystemExit(main())
