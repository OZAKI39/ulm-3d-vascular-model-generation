@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  py -3.14 -m venv .venv
  if errorlevel 1 goto fail
)
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto fail
echo Ready. Double-click start_chronos.cmd.
pause
exit /b 0
:fail
echo Setup failed. Check Python 3.14 and the internet connection.
pause
exit /b 1
