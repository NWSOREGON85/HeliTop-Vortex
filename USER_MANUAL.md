HeliTop Vortex v15.0 - User Manual (Plain Text Version)
============================================================

Version: 15.0
Date: March 2026
Developed by: Nathaniel & Grok (xAI)

1. INTRODUCTION
---------------
HeliTop Vortex is a fast, free, GPU-accelerated vortex-filament simulator 
with an intuitive graphical user interface. It uses the Biot-Savart law, 
adaptive regridding, and topology-aware reconnection to simulate vortex 
wakes for marine propellers, rocket plumes, aircraft wakes, circular pipes, 
and generic flows.

Key advantages:
- Real-time Live 3D Preview of vortex filaments
- Multiple engineering presets
- Cylinder, flat-wall confinement, and buoyancy/thermal plumes
- GPU acceleration (CuPy) with automatic CPU fallback
- Professional outputs: rotating 3D GIFs, PDF reports, dynamic CSV files, VTK files
- Built-in validation suite and one-click unit tests
- Clean, modern GUI with sliders and console log

2. INSTALLATION
---------------
Windows (Recommended - Single EXE):
- Download HeliTop_Vortex_v15.0.exe from the Releases page
- Double-click to run - no installation required

From Source (Advanced):
- git clone https://github.com/yourusername/heli-top-vortex.git
- cd heli-top-vortex
- pip install numpy matplotlib pillow vtk cupy-cuda12x   (or cupy-cuda13x)
- python heli_top_gui.py

3. MAIN INTERFACE GUIDE - EVERY CONTROL EXPLAINED
-------------------------------------------------

TOP BAR CONTROLS
----------------
• Preset (dropdown menu)
  Options: marine_propeller, rocket_plume, marine_propeller_high_skew, aircraft_wake, generic, circular_pipe
  What it does: Chooses the type of simulation. Each preset automatically loads realistic filament geometry and circulation values.
  Example: Choose "rocket_plume" to simulate a pulsating rocket exhaust.

• Confinement (dropdown menu)
  Options: none, cylinder, flat_wall
  What it does: Adds physical boundaries.
  - "none" → free space (default)
  - "cylinder" → simulates flow inside a pipe (uses image vortices)
  - "flat_wall" → simulates flow near a flat surface (uses image vortices)
  Example: Set Confinement = cylinder and preset = circular_pipe to model swirling flow inside a pipe.

• Buoyancy (checkbox)
  What it does: Enables thermal buoyancy forces (useful for rocket plumes or hot gas rises).
  When checked, the buoyancy strength slider (below) becomes active.

• GPU (checkbox)
  What it does: Turns on CuPy GPU acceleration. If your computer has an NVIDIA GPU with CUDA, the simulation runs much faster.
  Recommended: Leave checked unless you get errors.

• GIF (checkbox)
  What it does: Automatically saves rotating 3D GIF animations of every realization.
  Uncheck if you only want CSV/PDF and don't need the GIFs (saves disk space).

• Dark (checkbox)
  What it does: Switches the entire GUI to dark mode (black background, green console text). Click to toggle.

PROGRESS BAR (top of window)
----------------------------
Shows overall campaign progress (0-100%) across all realizations and steps. Updates live.

LEFT SIDE - MAIN PARAMETERS SLIDERS
-----------------------------------
• N_FIL (points per filament)
  Range: 64 to 512
  What it does: Controls resolution of each vortex filament. Higher = smoother but slower.
  Example: Set to 384 for high-quality propeller wake studies.

• Realizations (NUM_REALIZATIONS)
  Range: 1 to 10
  What it does: Runs multiple independent simulations for statistical results (Monte-Carlo style).
  Example: Set to 3 to get mean and standard deviation of thrust/efficiency.

• Steps
  Range: 50 to 300
  What it does: Number of time steps in each realization. More steps = longer physical time simulated.
  Example: Set to 150 for a full propeller revolution.

• Core Size
  Range: 0.01 to 0.2
  What it does: Base vortex core radius (regularization parameter). Prevents singularities.
  Example: Lower values (0.04) give sharper vortices.

• Pipe Radius
  Range: 0.5 to 5.0
  What it does: Only active for circular_pipe preset. Sets the pipe wall radius.
  Example: Set to 2.5 for a typical marine ducted propeller simulation.

• Swirl Layers
  Range: 1 to 8
  What it does: Only active for circular_pipe preset. Number of concentric vortex layers.
  Example: Set to 5 for realistic swirling pipe flow.

• Axial Strength
  Range: 0.5 to 2.0
  What it does: Strength of axial (streamwise) velocity component in pipe flows.
  Example: Set to 1.8 for strong forward flow in a pipe.

• Friction
  Range: 0.01 to 0.1
  What it does: Wall friction factor for pipe flows.
  Example: Set to 0.03 for realistic pipe wall drag.

• Buoyancy
  Range: 0.1 to 2.0
  What it does: Strength of thermal buoyancy force (only active when Buoyancy checkbox is checked).
  Example: Set to 0.8 for moderate plume rise in rocket exhaust.

RIGHT SIDE - THRUST CALIBRATION SLIDERS
---------------------------------------
• Circulation ×
  Range: 0.5 to 2.0
  What it does: Multiplies the strength of all vortices.
  Example: Set to 1.3 to increase thrust by ~30%.

• Core Size ×
  Range: 0.5 to 2.0
  What it does: Multiplies the base core size.
  Example: Set to 0.8 for tighter, more intense vortices.

• Viscosity ×
  Range: 0.5 to 2.0
  What it does: Multiplies the stochastic noise (viscous diffusion).
  Example: Set to 1.2 for slightly more turbulent diffusion.

• Thrust Scale ×
  Range: 0.5 to 2.0
  What it does: Scales all computed thrust and torque values.
  Example: Set to 1.15 to match experimental thrust data.

• Enstrophy Cap
  Range: 1e6 to 5e7
  What it does: Caps maximum enstrophy to prevent numerical blow-up.
  Example: Set to 5e6 for very long runs.

• Dynamic Amplitude
  Range: 0.0 to 0.2
  What it does: Adds sinusoidal oscillation to thrust (rocket_plume only).
  Example: Set to 0.1 for ±10% thrust pulsing.

RIGHT COLUMN BUTTONS
--------------------
• RUN
  Starts the full simulation campaign with current settings.

• Validate
  Runs a quick physics validation suite (vortex ring + leapfrogging rings) and saves a PDF report. Use this before long runs.

• Tests
  Runs all internal unit tests and shows results in the console. Confirms everything is working.

• Close
  Saves current settings to config.json and exits the program safely.

CENTER - LIVE 3D PREVIEW
------------------------
Shows real-time 3D vortex filaments updating automatically every 12 steps. No extra button needed.

BOTTOM - CONSOLE LOG
--------------------
Displays step-by-step progress, reconnection events, and any warnings (green text on black).

4. CLEAR EXAMPLES OF HOW TO USE THE APPLICATION
-----------------------------------------------

EXAMPLE 1: Marine Propeller Performance Study
---------------------------------------------
1. Select Preset = marine_propeller
2. Set N_FIL = 384, Realizations = 3, Steps = 150
3. Set Circulation × = 1.25, Thrust Scale × = 1.1
4. Leave Confinement = none, GPU checked, GIF checked
5. Click RUN
Result: You get Kt/Kq/η values in the CSV, efficiency curve in the PDF, and rotating 3D GIFs.

EXAMPLE 2: Pulsating Rocket Plume
---------------------------------
1. Select Preset = rocket_plume
2. Set Dynamic Amplitude = 0.12, Thrust Scale × = 1.15
3. Check Buoyancy checkbox and set Buoyancy slider = 0.8
4. Set Steps = 120, Realizations = 2
5. Click RUN
Result: Thrust history plot in PDF, pulsating plume visible in Live 3D Preview and GIF.

EXAMPLE 3: Confined Pipe Flow
-----------------------------
1. Select Preset = circular_pipe
2. Set Confinement = cylinder
3. Set Pipe Radius = 2.5, Swirl Layers = 5, Axial Strength = 1.8
4. Set Friction = 0.03
5. Click RUN
Result: Realistic swirling flow inside a pipe with wall effects, saved as VTK files for ParaView.

5. USING THE config.json FILE - STEP-BY-STEP
--------------------------------------------
1. Run the program once (it creates config.json automatically)
2. Close the program
3. Open config.json in Notepad
4. Change any value (e.g., "dynamic_amplitude": 0.15)
5. Save the file
6. Restart the program
7. Click RUN - your new settings are loaded automatically

6. TIPS & BEST PRACTICES
------------------------
- Start with small settings (N_FIL=128, Realizations=1) for quick tests
- Always click "Validate" first on a new preset
- Enable GPU for faster runs
- Use GIF checkbox only when you need animations

7. TROUBLESHOOTING
------------------
Problem: Live 3D Preview is blank → Make sure you are using v15.0
Problem: Out of memory → Reduce N_FIL or Realizations
Problem: CuPy warning → Reinstall the correct CuPy version for your CUDA

Thank you for using HeliTop Vortex!
Feedback and feature requests are always welcome.

End of Manual
