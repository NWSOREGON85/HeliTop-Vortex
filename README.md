# HeliTop Vortex v13.0

**A fast, free, GPU-accelerated vortex-filament simulator with an intuitive GUI**  
Perfect for marine propellers, rocket plumes, aircraft wakes, and general vortex-dynamics studies.

![HeliTop Vortex Screenshot](https://github.com/yourusername/helitop-vortex/screenshots/gui-main.png)  
*(Replace with your actual screenshot after uploading)*

Why HeliTop Vortex?

HeliTop Vortex is a **mid-fidelity Lagrangian vortex-filament tool** designed for engineers and researchers who need quick, physically realistic wake simulations **without** the complexity and cost of full CFD.

It combines:
- Stochastic noise (Navier–Stokes inspired viscosity)
- Adaptive regridding
- Topology-aware reconnection
- Real-time calibration sliders
- Beautiful 3D rotating GIF animations

All wrapped in a clean, single-file Windows EXE that runs instantly — **no installation required**.

Features

- **Modern Tkinter GUI** with live calibration sliders (circulation, core size, viscosity, thrust scale, dynamic amplitude, enstrophy cap)
- **GPU acceleration** via CuPy (automatic fallback to CPU/NumPy)
- **Stochastic Lagrangian vortex filaments** with Biot-Savart induction
- **Adaptive regridding** and **topology-based reconnection**
- **Built-in presets**:
  - Marine Propeller
  - Rocket Plume (with time-varying thrust)
  - Marine Propeller (High Skew)
  - Aircraft Wake
  - Generic (custom rings)
- **Rich outputs**:
  - 3D rotating GIF animations (per realization)
  - Professional PDF reports
  - VTK files for ParaView
  - CSV performance data
  - Efficiency curves (Kt/Kq/η) and thrust profiles
- **Config persistence** via `config.json`
- **Dark mode** support
- **Fully self-contained** single EXE (Windows)

Screenshots & Demo

*(Upload your best screenshots and GIFs here after creating the repo)*

- Main GUI with calibration panel
- Marine Propeller efficiency curve
- Rocket Plume thrust profile
- 3D rotating vortex animation (example GIF)

Installation

### Option 1: Windows EXE (Recommended)

1. Download the latest `HeliTop_Vortex_v13.0.exe` from the [Releases page](https://github.com/yourusername/heli-top-vortex/releases)
2. Run it — no installation needed

### Option 2: From Source

```bash
git clone https://github.com/yourusername/heli-top-vortex.git
cd heli-top-vortex
pip install numpy matplotlib pillow vtk cupy-cuda12x  # or cupy-cuda13x
python heli_top_gui.py
