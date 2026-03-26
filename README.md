# HeliTop Vortex v15.0

**Fast, free, GPU-accelerated vortex-filament simulator with intuitive GUI**

A mid-fidelity Lagrangian vortex method tool designed for quick wake-flow studies in marine propulsion, rocket plumes, aircraft wakes, and confined flows (pipes). Built for engineers who need results in seconds to minutes instead of hours.

**Version:** 15.0 (March 2026)  
**Developed by:** Nathaniel

![HeliTop Vortex Screenshot](https://github.com/NWSOREGON85/helitop-vortex/raw/main/screenshots/gui.png)

## ✨ Key Features

- **Live 3D Preview** – Real-time rotating vortex filaments during simulation
- **Multiple Engineering Presets** – Marine propeller, rocket plume, aircraft wake, circular pipe, generic
- **Confinement & Physics** – Cylinder, flat-wall, buoyancy/thermal plumes
- **GPU Acceleration** – CuPy (CUDA 12.x/13.x) with automatic CPU fallback
- **Professional Outputs** – 3D rotating GIFs, PDF reports, dynamic CSV performance data, VTK files for ParaView
- **Built-in Validation Suite** – Vortex ring + leapfrogging rings with conservation checks
- **Unit Tests** – One-click full system validation
- **Clean GUI** – Sliders, dark mode, progress bar, console log
- **JSON Configuration** – Save/load all settings

## Installation (Windows – Recommended)

1. Download the latest release: `HeliTop_Vortex_v15.0.exe`
2. Run the EXE – **no installation required**
3. (Optional) Place in any folder; it will create `plots/`, `reports/`, `vtk/`, and `config.json` automatically.

### From Source

```bash
git clone https://github.com/yourusername/heli-top-vortex.git
cd heli-top-vortex
pip install numpy matplotlib pillow vtk cupy-cuda12x  # or cupy-cuda13x
python heli_top_gui.py
