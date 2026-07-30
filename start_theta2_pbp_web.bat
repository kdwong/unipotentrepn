@echo off
setlocal
title Theta Atlas - BMSZ painted bipartitions
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
  echo The prepared Python environment was not found in:
  echo   %CD%\.venv
  echo.
  echo Recreate the environment before starting the webpage.
  pause
  exit /b 1
)

echo Starting Theta Atlas...
echo Keep this window open while using the webpage.
echo.
".venv\Scripts\python.exe" "theta2_pbp_web.py" %*

if errorlevel 1 (
  echo.
  echo The painted-bipartition webpage could not be started.
  echo Check the message above, then press any key to close this window.
  pause >nul
)

endlocal
