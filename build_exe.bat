@echo off
pip install pyinstaller numpy matplotlib --quiet
pyinstaller --onefile --windowed --hidden-import=tkinter --hidden-import=matplotlib.backends.backend_tkagg --name=HeliTop_Vortex heli_top_gui.py
echo.
echo ✅ EXE created in dist\HeliTop_Vortex.exe
echo Double-click it now — it should open the GUI.
pause
