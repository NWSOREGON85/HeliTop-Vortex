@echo off
pip install pyinstaller --quiet
pyinstaller --onefile --windowed --name=HeliTop_Vortex heli_top_gui.py
echo ✅ EXE created in dist/ folder!
pause
