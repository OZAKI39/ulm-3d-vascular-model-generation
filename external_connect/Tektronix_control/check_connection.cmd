@echo off
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
  ".venv\Scripts\python.exe" -X utf8 afg_control.py discover
) else (
  python -X utf8 afg_control.py discover
)
echo Press any key to close this window...
pause >nul
