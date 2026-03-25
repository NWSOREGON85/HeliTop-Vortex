@echo off
echo ================================================
echo Building HeliTop Vortex v12.4 EXE - Final Release
echo ================================================

mkdir plots 2>nul
mkdir reports 2>nul
mkdir vtk 2>nul

pyinstaller --onefile --console --name "HeliTop_Vortex_v12.4" ^
  --add-data "plots;plots" ^
  --add-data "reports;reports" ^
  --add-data "vtk;vtk" ^
  --hidden-import tkinter ^
  --hidden-import matplotlib ^
  --hidden-import PIL ^
  --hidden-import numpy ^
  --hidden-import cupy ^
  heli_top_gui.py

echo.
echo ================================================
echo Build finished!
echo EXE is in the "dist" folder.
echo ================================================
pause
