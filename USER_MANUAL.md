HeliTop Vortex v13.0 — User Manual

A fast, free, GPU-accelerated vortex-filament simulator with an intuitive GUI

Version: 13.0 Final Release
Date: March 2026
Developed by: Nathaniel & Grok (xAI)

Table of Contents
1. Introduction
2. Installation
3. Quick Start
4. Main Interface Guide
5. Presets
6. Calibration Panel
7. Using the config.json File
8. Understanding the Outputs
9. Tips & Best Practices
10. Troubleshooting
11. Roadmap & Contributing

Introduction
HeliTop Vortex is a mid-fidelity Lagrangian vortex-filament simulator designed for quick wake-flow studies in marine propulsion, rocket plumes, and aircraft design. It uses stochastic noise, adaptive regridding, and topology-aware reconnection to produce physically realistic results in seconds to minutes instead of hours.

Key advantages:
- Modern, easy-to-use GUI (no command line required)
- GPU acceleration (CuPy) with automatic CPU fallback
- Real-time calibration sliders
- Professional outputs: 3D rotating GIFs, PDF reports, VTK files, CSV data

Installation

Windows (Recommended — Single EXE)
1. Download HeliTop_Vortex_v13.0.exe from the Releases page
2. Run the EXE — no installation needed

From Source (Advanced)
git clone https://github.com/yourusername/heli-top-vortex.git
cd heli-top-vortex
pip install numpy matplotlib pillow vtk cupy-cuda12x   # or cupy-cuda13x
python heli_top_gui.py

Quick Start
1. Launch the program
2. Select a preset (e.g., marine_propeller)
3. Adjust any sliders if desired
4. Click RUN SIMULATION
5. Watch the live progress bar and console
6. When finished, check the folders:
   - reports/ → PDF report
   - plots/ → efficiency curve + 3D GIFs
   - vtk/ → files for ParaView
   - propeller_performance.csv

Main Interface Guide
- Preset dropdown — choose simulation type
- Checkboxes:
  - Use GPU Acceleration (CuPy)
  - Save Vortex Animation as GIF
  - Dark Mode
- Parameters section — N_FIL, Realizations, Steps, Core Size
- Thrust Calibration panel — detailed below
- Progress bar & live console — real-time feedback
- RUN SIMULATION button
- Close Application button

Presets
marine_propeller        Ship propellers                   Kt/Kq/η curves
rocket_plume            Rocket exhaust plumes             Thrust time-history
marine_propeller_high_skew   High-skew propellers      Higher efficiency cases
aircraft_wake           Wing-tip vortices                 Vortex pair evolution
generic                 Custom ring vortices              General testing

Calibration Panel (Rocket Plume & Propellers)
All multipliers are relative (1.0 = default):

- Circulation × — vortex strength (higher = stronger vortices)
- Core Size × — smoothing radius (higher = more stable, less sharp)
- Viscosity × — stochastic noise level (higher = more diffusion)
- Thrust Scale × — overall force scaling
- Enstrophy Cap — maximum allowed enstrophy (for stability)
- Dynamic Amplitude — only for rocket_plume; adds sinusoidal thrust variation (0.0 = steady, 0.1–0.2 = realistic pulsing)

Example: For a stronger rocket plume, set Circulation × = 1.4 and Dynamic Amplitude = 0.12.

Using the config.json File
The program automatically creates and updates config.json every time you run a simulation.

How to edit it manually (advanced users)
1. Close the program
2. Open config.json in any text editor (Notepad, VS Code, etc.)
3. Example content:

{
  "preset": "rocket_plume",
  "N_FIL": 256,
  "NUM_REALIZATIONS": 3,
  "steps": 150,
  "core_base": 0.08,
  "use_gpu": true,
  "save_gif": true,
  "dark_mode": true,
  "circulation_multiplier": 1.35,
  "core_multiplier": 1.1,
  "viscosity_multiplier": 0.9,
  "thrust_scale_multiplier": 1.2,
  "enstrophy_cap": 15000000,
  "dynamic_amplitude": 0.15
}

4. Save the file
5. Restart the program — your settings will be loaded automatically

Tip: Use this file to create “preset configurations” (e.g., one for Starship plume studies, one for propeller optimization).

Understanding the Outputs
- reports/HeliTop_*.pdf — Professional summary with all key metrics
- plots/efficiency_curve.png — Plot of efficiency or thrust vs. time/step
- plots/vortex_animation_r1_3D.gif — Rotating 3D vortex visualization (one per realization)
- vtk/plume_step_XXXX.vtk — Load in ParaView for interactive 3D viewing
- propeller_performance.csv — Raw data table (easy to open in Excel)

Tips & Best Practices
- Start with 2 realizations and 100 steps for quick tests
- Increase N_FIL to 512 for higher fidelity (slower)
- For rocket plumes, enable Dynamic Amplitude between 0.08–0.18
- Save GIFs only when you need animations (they take extra time)
- Use ParaView to open VTK files for beautiful 3D animations

Troubleshooting
Problem: Program closes immediately
Solution: Run build_exe.bat again or run heli_top_gui.py from command prompt to see errors.

Problem: Console is gray instead of black
Solution: Make sure you are using the latest v13.0 file.

Problem: Out of memory
Solution: Reduce N_FIL to 128 or turn off GPU.

Problem: No progress bar movement
Solution: Rebuild the EXE with the latest heli_top_gui.py

Roadmap & Contributing
Future features planned:
- Import experimental wake data
- Export to ANSYS/Star-CCM+
- Batch mode for multiple runs
- Linux & macOS support

Want to contribute?
Fork the repo, make changes, and open a Pull Request. All contributions are welcome!

Thank you for using HeliTop Vortex!
