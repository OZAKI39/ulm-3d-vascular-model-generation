"""Open records on WSL without relying on an installed Linux desktop launcher."""
import os
from pathlib import Path
import shutil
import subprocess


def is_wsl():
    return bool(os.environ.get('WSL_INTEROP') or os.environ.get('WSL_DISTRO_NAME'))


def open_directory(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=True)
    errors = []
    if is_wsl():
        try:
            converter = shutil.which('wslpath') or '/usr/bin/wslpath'
            explorer = shutil.which('explorer.exe') or '/mnt/c/Windows/explorer.exe'
            windows_path = subprocess.run([converter, '-w', str(path)], check=True,
                capture_output=True, text=True, timeout=3).stdout.strip()
            if not windows_path:
                raise OSError('wslpath returned an empty path')
            # Pass the UNC path as one literal argument; never use a shell.
            subprocess.Popen([explorer, windows_path], stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
            return True, 'Opened records in Windows File Explorer.'
        except (OSError, subprocess.SubprocessError) as exc:
            errors.append(str(exc))
    from PySide6.QtCore import QUrl
    from PySide6.QtGui import QDesktopServices
    if QDesktopServices.openUrl(QUrl.fromLocalFile(str(path))):
        return True, 'Opened records folder.'
    return False, '\n'.join(errors + ['No file manager could be launched.', str(path)])
