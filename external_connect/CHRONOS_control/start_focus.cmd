@echo off
cd /d "%~dp0"
".venv\Scripts\python.exe" "%~dp0focus_assist.py"
if errorlevel 1 pause
