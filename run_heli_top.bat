@echo off
title HeliTop Vortex v8.0
echo Installing dependencies...
pip install numpy matplotlib --quiet
echo Launching GUI...
python heli_top_gui.py
pause
