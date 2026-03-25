# Quick Start Guide – HeliTop Vortex v13.0

**Get your first vortex simulation running in under 60 seconds.**

### 1. Download the program
- Go to the **[Releases](https://github.com/yourusername/heli-top-vortex/releases)** page
- Download the file: `HeliTop_Vortex_v13.0.exe`

### 2. Run the program
- Double-click the downloaded `.exe` file  
  (No installation required — it runs immediately)

### 3. Run your first simulation
1. In the **Preset** dropdown, choose **marine_propeller** (best for first run)
2. Leave all sliders at their default values
3. Click the big blue **RUN SIMULATION** button
4. Wait for the progress bar to reach 100% (usually 30–90 seconds)

### 4. Find your results
After the simulation finishes, the following folders and files appear in the same folder as the EXE:

- `reports/` → Professional PDF report with all key numbers
- `plots/` → Efficiency curve image + 3D rotating GIF animations
- `vtk/` → Files you can open in ParaView
- `propeller_performance.csv` → Raw data table

**You’re done!** You have successfully run a complete vortex simulation.

---

### Next steps (optional)
- Try the **rocket_plume** preset and move the **Dynamic Amplitude** slider
- Increase **N_FIL** to 512 for higher detail
- Check the **Dark Mode** box for a nicer look
- Edit `config.json` to save your favorite settings

Full documentation:  
- [README.md](README.md)  
- [USER_MANUAL.md](USER_MANUAL.md)

Enjoy exploring vortex dynamics!
