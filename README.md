# HeliTop Vortex v8.0

**The best free professional vortex filament simulator for engineers**

Fast, accurate, and easy-to-use tool for marine propellers, aircraft wakes, wind turbines, jet plumes, and more.

Built for real engineering work: propeller design, wake turbulence studies, vortex breakdown analysis, and rapid concept testing.

![HeliTop GUI Screenshot](https://via.placeholder.com/800x400/0066cc/ffffff?text=HeliTop+Vortex+GUI)  
*(Replace with your actual screenshot after first run)*

## Features
- Beautiful, modern GUI (no coding required)
- Marine propeller preset with automatic **Kt / Kq** and **efficiency curves**
- Aircraft wake preset with induced drag proxies
- Dynamic topology-triggered reconnection (prevents singularities)
- One-click PDF engineering reports
- ParaView-ready VTK export + live reconnection movie (GIF)
- Efficiency vs Advance Ratio (J) curves
- One-click .exe version (no Python install needed)

## Quick Start (Windows – 30 seconds)

1. Download the latest release ZIP
2. Extract the folder
3. Double-click **`run_heli_top.bat`**
4. Choose your preset (e.g., "marine_propeller")
5. Click **RUN SIMULATION**

You’ll instantly get:
- `propeller_performance.csv`
- `efficiency_curve.png`
- `reports/HeliTop_Propeller_Report.pdf`
- VTK files for ParaView
- Reconnection animation GIF

## One-Click Executable
Double-click **`build_exe.bat`** to create a single `HeliTop_Vortex.exe` (no Python required).

## Screenshots

*(Add 3–4 screenshots here after your first run: GUI, efficiency curve, ParaView visualization, reconnection movie)*

## Supported Presets
- `marine_propeller` (default – full Kt/Kq + efficiency)
- `aircraft_wake`
- `generic` (custom flows)

## License
MIT License – free for commercial and academic use.

## Acknowledgments
Developed as a collaborative engineering tool. Inspired by real-world needs in marine propulsion and aerospace wake analysis.

---

**Star this repo** if you find it useful!  
Questions or want a new preset (wind turbine, rocket plume, etc.)? Open an Issue or PR.
