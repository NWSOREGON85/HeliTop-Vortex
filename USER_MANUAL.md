# HeliTop Vortex v18.0 — User Manual

A fast, free, GPU-accelerated vortex-filament simulator with an intuitive GUI

Version: 18.0 Final Release  
Date: March 2026  
Developed by: Nathaniel

## Table of Contents
1. Introduction  
2. Installation  
3. Quick Start  
4. Main Interface Guide  
5. Presets  
6. Calibration Panel (Sliders & Controls)  
7. Using the config.json File  
8. Understanding the Outputs  
9. Tips & Best Practices  
10. Troubleshooting  
11. Roadmap & Contributing  

## Introduction
HeliTop Vortex is a mid-fidelity Lagrangian vortex-filament simulator designed for quick wake-flow studies in marine propulsion, rocket plumes, aircraft wakes, HVAC ducts, and circular pipe flows. It uses adaptive vortex filaments, stochastic diffusion, reconnection, and image vortices to produce physically realistic results in seconds to minutes on a normal laptop.

**Best academic uses**  
In academia the tool is excellent for research in vortex dynamics, topological fluid mechanics, and Lagrangian methods. It is particularly useful for exploring depletion mechanisms and linking growth in Navier-Stokes equations, for validating new turbulence models, and for educational purposes in computational fluid dynamics courses. Students and researchers can quickly test hypotheses about vortex reconnection, wall shear, pressure drop, and efficiency without needing expensive supercomputing resources.

**Best professional uses**  
Professionally, it is highly valued by engineers working on HVAC system design (where accurate wall shear stress and pressure drop predictions are critical for fan sizing and noise reduction), marine propulsion (propeller and ducted thruster optimization), aerospace (aircraft wake modeling for airport safety and formation flight), and rocket propulsion (plume analysis for launch vehicles). The fast turnaround time and easy-to-use GUI make it ideal for rapid design iteration, concept screening, and preliminary studies before moving to full CFD in ANSYS or Star-CCM+.

Key advantages:
- Modern easy-to-use GUI with no command line required
- GPU acceleration with automatic CPU fallback
- Real-time live 3D preview plus rotating GIF animations
- Professional outputs: PDF reports, efficiency curves, wall-shear plots, VTU files for ParaView, and CSV performance data
- Built-in validation suite and unit tests

## Installation
For Windows (recommended):  
Download HeliTop_Vortex_v18.0.exe from the releases page and run it. No installation is needed.

For advanced users:  
Clone the repository, install numpy matplotlib pillow vtk and the correct cupy-cuda package for your GPU, then run python heli_top_gui.py.

## Quick Start
1. Launch the program.  
2. Choose a preset from the drop-down menu (for example circular_pipe or hvac_pipe).  
3. Leave the sliders at default values unless you want to experiment.  
4. Click the large RUN button.  
5. Watch the live 3D preview window and the console output.  
6. When the simulation finishes, check the folders:  
   - reports folder contains the full PDF report  
   - plots folder contains the efficiency curve PNG and rotating GIF animations  
   - vtk folder contains velocity field files you can open in ParaView

## Main Interface Guide
At the top you will see:  
- Preset drop-down menu: Choose the type of simulation (marine_propeller, rocket_plume, aircraft_wake, circular_pipe, hvac_pipe, etc.).  
- Confinement drop-down menu: none, cylinder, or flat_wall. Use cylinder for pipe or duct flows.  
- Checkboxes for GPU, GIF output, Dark mode, and Buoyancy.

On the left side are the main parameter sliders and controls:  
- N_FIL slider: Number of points per filament (64 to 512). Higher values give smoother filaments but take longer to run.  
- Realizations slider: Number of independent runs averaged together (1 to 10). Use 4 or more when you need reliable statistics.  
- Steps slider: Number of time steps in the simulation (50 to 300). 100 is the usual good balance.  
- Core Size slider: Vortex core radius (0.01 to 0.2). Our sweeps showed 0.10 gives the best efficiency.  
- Pipe Radius slider: Only active for pipe presets. Sets the radius of the duct.  
- Swirl Layers slider: Number of concentric rings of vorticity (1 to 8). Our sweeps showed 2 layers gives the highest efficiency.  
- Axial Strength slider: Strength of the main flow through the pipe (0.5 to 2.0). Higher values increase velocity and pressure drop.  
- Friction slider: Wall friction factor (0.01 to 0.1). Used only for HVAC pressure-drop calculations.  
- Buoyancy Strength slider: Strength of buoyancy force when buoyancy is enabled.  
- RPM slider: Rotational speed for propeller presets (100 to 3000). Only affects marine_propeller presets.

On the right side is the Calibration panel with extra multipliers:  
- Circulation multiplier slider: Scales overall vortex strength (0.5 to 2.0).  
- Core multiplier slider: Further scales the core radius (0.5 to 2.0).  
- Viscosity multiplier slider: Scales diffusion (0.5 to 2.0). Default 0.5 is optimal from our sweeps.  
- Thrust Scale multiplier slider: Scales the final thrust or force output.  
- Enstrophy Cap slider: Maximum allowed enstrophy to prevent numerical blow-up.  
- Dynamic Amplitude slider: Adds sinusoidal variation over time (useful for rocket plumes).

Below the sliders are four big buttons:  
- RUN: Starts the simulation.  
- Validate: Runs the full physics validation suite.  
- Tests: Runs quick unit tests on all features.  
- Close: Saves settings and exits cleanly.

## Presets
circular_pipe  
Best for ducted fans, marine thrusters, and wind-tunnel tests.  
Recommended settings: Core Size 0.10, Swirl Layers 2.  
Example: Design a high-efficiency ducted propeller for a drone. Set Swirl Layers to 2 and run 4 realizations. You will get efficiency around 0.242 with the lowest torque.

hvac_pipe  
Best for HVAC duct design and ventilation systems.  
Key outputs: wall shear stress distribution plot, pressure drop ΔP, and Reynolds number.  
Example 1 (low-velocity office duct): Pipe Radius 0.3 m, Axial Strength 1.0, Swirl Layers 3. Check that wall shear stays below 0.1 Pa and ΔP is under 1 Pa.  
Example 2 (high-velocity industrial exhaust): Pipe Radius 0.8 m, Axial Strength 1.8, Swirl Layers 2. Use the shear plot to decide whether you need a stronger fan.  
Example 3 (quiet HVAC system): Increase Core Size to 0.12 and Swirl Layers to 4. This gives smoother flow and lower peak shear for reduced noise.

marine_propeller  
Best for ship propellers and underwater vehicles.  
Example: Test a 4-blade propeller at 1200 RPM. Move the RPM slider live and watch how efficiency changes with advance ratio.

rocket_plume  
Best for rocket exhaust simulation.  
Example: Simulate a Starship-like plume by setting Dynamic Amplitude to 0.1. The report will show peak thrust and thrust oscillation over time.

## Calibration Panel – Detailed Slider Examples
**Core Size slider (0.01–0.2)**  
Controls the vortex core radius. Smaller values produce sharper, more intense vortices and usually higher efficiency, but they can lead to numerical instability if too small. Larger values add more diffusion, making the flow smoother and more stable but slightly lowering efficiency.  
HVAC example: If your duct shows excessively high wall shear stress, increase Core Size to 0.12 to dampen the vortices near the wall and reduce noise and energy losses.

**Swirl Layers slider (1–8)**  
Determines how many concentric rings of vorticity are placed inside the pipe. Fewer layers concentrate the swirl into a tighter core, giving higher efficiency and lower torque. More layers spread the vorticity outward, producing a smoother velocity profile but slightly higher torque and wall interaction.  
HVAC example: For a long straight ventilation duct that needs uniform flow, try 3 or 4 layers. For maximum efficiency in a short fan section, use only 2 layers.

**Axial Strength slider (0.5–2.0)**  
Controls the strength of the main axial flow through the pipe. Higher values increase bulk velocity and therefore pressure drop and Reynolds number.  
HVAC example: Set to 1.0 for normal office air-conditioning ducts. Set to 1.8 when designing a high-speed industrial exhaust system that needs strong ventilation.

**Friction slider (0.01–0.1)**  
Sets the wall friction factor used in pressure-drop calculations for HVAC and pipe presets.  
HVAC example: Use 0.02 for smooth metal ducts. Increase to 0.05–0.08 when modeling rough concrete or lined ducts.

**RPM slider (100–3000)**  
Sets rotational speed for propeller presets only. Higher RPM increases thrust and torque.  
Example: Set 1800 RPM to test a high-speed marine thruster and watch how efficiency changes with advance ratio.

**Viscosity multiplier slider (0.5–2.0)**  
Scales the amount of viscous diffusion. Our sweeps showed 0.5 gives the best balance of realism and efficiency.  
HVAC example: Increase to 1.0 or higher if you want a more diffusive, smoother flow field (useful for preliminary studies).

**Dynamic Amplitude slider (0.0–0.2)**  
Adds sinusoidal time-varying thrust (primarily for rocket plume presets).  
Example: Set to 0.1 to simulate realistic plume oscillation during launch.

## Using the config.json File
The program automatically saves and loads a file called config.json in the same folder as the executable.

How to use it:  
1. Run the program once. It creates config.json.  
2. Close the program.  
3. Open config.json in Notepad or any text editor.  
4. Change any value you want (for example change "num_swirl_layers" to 2 or "core_base" to 0.10).  
5. Save the file.  
6. Restart the program. Your new settings are automatically loaded.

Useful HVAC example inside config.json:  
Change the following lines:  
"preset": "hvac_pipe",  
"pipe_radius": 0.4,  
"num_swirl_layers": 3,  
"axial_strength": 1.5,  
"core_base": 0.10,  
"viscosity_multiplier": 0.5  

This gives a realistic HVAC duct with moderate velocity and the best efficiency we found.

## Understanding the Outputs
- PDF report in the reports folder: Contains all settings, mean Kt/Kq/efficiency (or thrust for plumes), wall shear stress, pressure drop, and Reynolds number.  
- efficiency_curve.png in the plots folder: Shows efficiency vs advance ratio for propellers, or wall shear distribution for HVAC.  
- Rotating GIF files in the plots folder: 3D animation of the vortex filaments.  
- VTU files in the vtk folder: 3D velocity fields you can open directly in ParaView for detailed inspection.  
- CSV file: Full data table for every realization (useful for spreadsheets).

## Tips & Best Practices
- Always start with the default settings. They are already optimized from our core-size and swirl-layer sweeps.  
- For HVAC design always look at the wall shear plot and the pressure drop value in the PDF.  
- Use 4 realizations when you need reliable average numbers.  
- Enable VTU export when you want to visualize the full 3D velocity field in ParaView.  
- If efficiency looks too low, try reducing swirl layers or slightly increasing core size.

## Troubleshooting
- No live preview appears: Make sure matplotlib is installed.  
- GPU checkbox does nothing: Install the correct cupy-cuda package for your graphics card.  
- Progress bar stays at 0 percent: The simulation is still running. Watch the console for step messages.  
- GIF creation fails: Close any other matplotlib windows before clicking RUN.

## Roadmap & Contributing
Current version (v18.0) already includes:  
- Optimal defaults from core-size & swirl-layer sweeps  
- Full HVAC pressure-drop & shear-stress analysis  
- VTU export for ParaView/Star-CCM+

Future plans: rotating bodies, wall boundary layers, ANSYS export.

Contributions are welcome! Fork the repo and send pull requests.

Enjoy designing better flows!  
Nathaniel
