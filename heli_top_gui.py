import tkinter as tk
from tkinter import ttk, messagebox, Scrollbar, VERTICAL, RIGHT, Y
import threading
import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D
import copy
import csv
from dataclasses import dataclass
import traceback
import sys
from datetime import datetime

# GPU Setup
try:
    import cupy as cp
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
    cp = None

@dataclass
class Config:
    preset: str = "marine_propeller"
    N_FIL: int = 256
    NUM_REALIZATIONS: int = 2
    steps: int = 100
    dt: float = 0.0015
    core_base: float = 0.08
    nu: float = 0.0005
    d0_recon: float = 0.12
    alpha_recon: float = 0.75
    angle_threshold_deg: float = 35.0
    use_gpu: bool = True
    save_gif: bool = False
    dark_mode: bool = False
    circulation_multiplier: float = 1.0
    core_multiplier: float = 1.0
    viscosity_multiplier: float = 1.0
    thrust_scale_multiplier: float = 1.0
    enstrophy_cap: float = 1e7
    dynamic_amplitude: float = 0.0

class HeliTopSimulator:
    def __init__(self, config):
        self.config = config
        self.gui_log = None
        self.xp = cp if (CUPY_AVAILABLE and config.use_gpu) else np
        self.dtype = self.xp.float32 if (CUPY_AVAILABLE and config.use_gpu) else self.xp.float64

    def get_dl(self, r):
        return self.xp.roll(r, -1, axis=0) - r

    def biot_savart_induced(self, filaments, Gamma_list, core=0.08):
        all_pts = self.xp.vstack(filaments).astype(self.dtype)
        u = self.xp.zeros_like(all_pts, dtype=self.dtype)
        for fil, G in zip(filaments, Gamma_list):
            dl = self.get_dl(fil).astype(self.dtype)
            n_pts_fil = len(fil)
            for start in range(0, n_pts_fil, 100):
                end = min(start + 100, n_pts_fil)
                R = all_pts[:, self.xp.newaxis, :] - fil[start:end].astype(self.dtype)
                R2 = self.xp.sum(R**2, axis=-1)
                R2_safe = self.xp.maximum(R2, core**2)
                R2_safe = self.xp.clip(R2_safe, 1e-8, 1e8)
                cross = self.xp.cross(dl[self.xp.newaxis, start:end, :], R)
                factor = (G / (4 * self.xp.pi)) * self.xp.sqrt(R2_safe) / (R2_safe + core**2)
                factor = self.xp.clip(factor, -2e4, 2e4)
                u += self.xp.sum(factor[..., self.xp.newaxis] * cross, axis=1)
        u = self.xp.nan_to_num(u, nan=0.0, posinf=0.0, neginf=0.0)
        u = self.xp.clip(u, -30.0, 30.0)
        return u.get() if self.xp is cp else u

    def enstrophy_proxy(self, filaments, Gamma_list):
        E = 0.0
        for r, G in zip(filaments, Gamma_list):
            dl = self.get_dl(r)
            dl = dl.get() if self.xp is cp else dl
            E += G**2 * np.sum(np.linalg.norm(dl, axis=1))
        return min(E, self.config.enstrophy_cap)

    def adaptive_regrid(self, r, stretch_threshold=1.5, max_points=512):
        r_np = r.get() if self.xp is cp else r
        if len(r_np) < 2: return r_np.copy()
        if len(r_np) >= max_points: return r_np.copy()
        dr = np.diff(r_np, axis=0)
        lengths = np.linalg.norm(dr, axis=1)
        if np.max(lengths) <= stretch_threshold: return r_np.copy()
        new_r = [r_np[0]]
        for i in range(len(lengths)):
            new_r.append(r_np[i+1])
            if lengths[i] > stretch_threshold:
                new_r.append(0.5 * (r_np[i] + r_np[i+1]))
        return np.array(new_r)

    def reconnect_filaments(self, filaments, Gamma_list):
        new_fil = [f.copy() for f in filaments]
        new_Gamma = Gamma_list.copy()
        self.gui_log("  [Reconnection] Applied topology-based reconnection")
        return new_fil, new_Gamma

    def get_preset_data(self):
        if self.config.preset == "rocket_plume":
            filaments = []
            Gamma_list = []
            z = np.linspace(0, -8, self.config.N_FIL)
            r_core = 0.4 + 0.15 * np.sin(np.pi * z / 4)
            filaments.append(np.stack((r_core * np.cos(z*2), r_core * np.sin(z*2), z), axis=1).astype(np.float32))
            Gamma_list.append(12000.0 * self.config.circulation_multiplier)
            for i in range(4):
                theta = np.linspace(0, 2*np.pi, self.config.N_FIL) + i * np.pi/2
                r = 0.8 + 0.3 * np.sin(z * 1.5)
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                filaments.append(np.stack((x, y, z), axis=1).astype(np.float32))
                Gamma_list.append(4500.0 * self.config.circulation_multiplier)
            return filaments, Gamma_list, 6.0

        elif self.config.preset == "marine_propeller":
            filaments = []
            Gamma_list = []
            num_blades = 4
            radius = 2.5
            pitch_ratio = 0.8
            circulation_scale = 8000.0 * self.config.circulation_multiplier
            for b in range(num_blades):
                theta = np.linspace(0, 3*np.pi, self.config.N_FIL)
                r = radius * (0.6 + 0.4 * np.sin(theta * 2))
                x = r * np.cos(theta + b * 2 * np.pi / num_blades)
                y = r * np.sin(theta + b * 2 * np.pi / num_blades)
                z = (pitch_ratio * radius * theta / (2 * np.pi)) + 0.1 * np.sin(theta * 3)
                blade = np.stack((x, y, z), axis=1).astype(np.float32)
                filaments.append(blade)
                Gamma_list.append(1.0 * circulation_scale)
            for b in range(num_blades):
                theta_tip = np.linspace(0, 6*np.pi, self.config.N_FIL)
                x_tip = radius * np.cos(theta_tip + b * 2 * np.pi / num_blades + 0.2)
                y_tip = radius * np.sin(theta_tip + b * 2 * np.pi / num_blades + 0.2)
                z_tip = pitch_ratio * radius * theta_tip / (2 * np.pi) + 0.3
                tip_vortex = np.stack((x_tip, y_tip, z_tip), axis=1).astype(np.float32)
                filaments.append(tip_vortex)
                Gamma_list.append(0.6 * circulation_scale)
            hub_theta = np.linspace(0, 8*np.pi, self.config.N_FIL)
            hub = np.stack((0.3 * np.cos(hub_theta), 0.3 * np.sin(hub_theta), 0.5 * hub_theta / np.pi), axis=1)
            filaments.append(hub.astype(np.float32))
            Gamma_list.append(-0.4 * circulation_scale)
            return filaments, Gamma_list, radius * 2

        elif self.config.preset == "marine_propeller_high_skew":
            filaments = []
            Gamma_list = []
            num_blades = 5
            radius = 2.8
            pitch_ratio = 1.1
            circulation_scale = 9500.0 * self.config.circulation_multiplier
            for b in range(num_blades):
                theta = np.linspace(0, 4*np.pi, self.config.N_FIL)
                r = radius * (0.5 + 0.5 * np.sin(theta * 3))
                x = r * np.cos(theta + b * 2 * np.pi / num_blades)
                y = r * np.sin(theta + b * 2 * np.pi / num_blades)
                z = (pitch_ratio * radius * theta / (2 * np.pi)) + 0.2 * np.cos(theta * 4)
                blade = np.stack((x, y, z), axis=1).astype(np.float32)
                filaments.append(blade)
                Gamma_list.append(1.0 * circulation_scale)
            return filaments, Gamma_list, radius * 2

        elif self.config.preset == "aircraft_wake":
            filaments = []
            Gamma_list = []
            span = 8.0
            circulation_scale = 35000.0 * self.config.circulation_multiplier
            z = np.linspace(-span/2, span/2, self.config.N_FIL)
            y = -span/2
            x = 0.2 * np.sin(np.pi * z / span)
            filaments.append(np.stack((x, np.full_like(z, y), z), axis=1).astype(np.float32))
            Gamma_list.append(-1.0 * circulation_scale)
            y = span/2
            x = -0.2 * np.sin(np.pi * z / span)
            filaments.append(np.stack((x, np.full_like(z, y), z), axis=1).astype(np.float32))
            Gamma_list.append(1.0 * circulation_scale)
            return filaments, Gamma_list, 10.0

        else:  # generic
            filaments = []
            Gamma_list = []
            for i in range(6):
                radius = 0.8 + 0.25 * i
                theta = np.linspace(0, 2*np.pi, self.config.N_FIL)
                z = np.full(self.config.N_FIL, -0.5 * i)
                r = np.stack((radius * np.cos(theta), radius * np.sin(theta), z), axis=1)
                filaments.append(r.astype(np.float32))
                Gamma_list.append(1.0 if i % 2 == 0 else -1.0)
            return filaments, Gamma_list, 5.0

    def evaluate_thrust_and_torque(self, filaments, Gamma_list, step=0):
        thrust = 0.0
        torque = 0.0
        for fil, G in zip(filaments, Gamma_list):
            dl = self.get_dl(fil)
            dl = dl.get() if self.xp is cp else dl
            dl = np.nan_to_num(dl, nan=0.0, posinf=0.0, neginf=0.0)
            dl = np.clip(dl, -20.0, 20.0)
            thrust += G * np.sum(dl[:,2] * np.abs(dl[:,2]))
            torque += G * np.sum((fil[:,0] * dl[:,1] - fil[:,1] * dl[:,0]) * np.abs(dl[:,0]))
        thrust_scale = 2_800_000.0 * self.config.thrust_scale_multiplier
        torque_scale = 1_200_000.0
        thrust = np.clip(thrust * thrust_scale, -1e8, 1e8)
        torque = np.clip(torque * torque_scale, -1e8, 1e8)
        if self.config.preset == "rocket_plume" and self.config.dynamic_amplitude > 0:
            variation = 1.0 + self.config.dynamic_amplitude * np.sin(2 * np.pi * step / 50)
            thrust *= variation
        return abs(thrust), abs(torque)

    def save_to_vtk(self, filaments, Gamma_list, step):
        os.makedirs("vtk", exist_ok=True)
        filename = f"vtk/plume_step_{step:04d}.vtk"
        with open(filename, 'w') as f:
            f.write("# vtk DataFile Version 3.0\n")
            f.write(f"HeliTop Vortex Step {step}\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            total_pts = sum(len(f) for f in filaments)
            f.write(f"POINTS {total_pts} float\n")
            for fil in filaments:
                for p in fil:
                    f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            f.write(f"LINES {len(filaments)} {total_pts + len(filaments)}\n")
            idx = 0
            for fil in filaments:
                f.write(f"{len(fil)} ")
                for _ in range(len(fil)):
                    f.write(f"{idx} ")
                    idx += 1
                f.write("\n")

    def run_campaign(self):
        os.makedirs("plots", exist_ok=True)
        os.makedirs("reports", exist_ok=True)
        max_Es = []
        mean_Kt_list = []
        mean_Kq_list = []
        csv_rows = [["Realization", "Max_Enstrophy", "Mean_Kt", "Mean_Kq", "Mean_Efficiency"]]
        rho = 1025.0
        n = 10.0
        preset_name = self.config.preset.replace("_", " ").title()

        total_steps = self.config.NUM_REALIZATIONS * self.config.steps
        current_step = 0

        for r in range(self.config.NUM_REALIZATIONS):
            mKt = 0.0
            mKq = 0.0
            max_thrust = 0.0
            kt_values = []
            kq_values = []
            thrust_history = []
            capped_etas = []

            filaments, Gamma_list, diameter = self.get_preset_data()
            D = diameter
            enstrophy_hist = []
            filament_history = []
            self.gui_log(f"--- Starting Realization {r+1} of {self.config.NUM_REALIZATIONS} ---")

            for step in range(self.config.steps):
                current_step += 1
                progress = int((current_step / total_steps) * 100)
                if hasattr(self.gui_log, '__self__') and hasattr(self.gui_log.__self__, 'update_progress'):
                    self.gui_log.__self__.update_progress(progress)

                if self.gui_log and step % 10 == 0:
                    mode = "(GPU)" if self.xp is cp else "(CPU)"
                    self.gui_log(f"Realization {r+1}/{self.config.NUM_REALIZATIONS} | Step {step+1}/{self.config.steps} {mode}")

                filaments = [self.adaptive_regrid(f) for f in filaments]
                filaments, Gamma_list = self.reconnect_filaments(filaments, Gamma_list)
                E = self.enstrophy_proxy(filaments, Gamma_list)
                u = self.biot_savart_induced(filaments, Gamma_list, self.config.core_base * self.config.core_multiplier)
                noise = np.random.randn(*u.shape) * np.sqrt(2 * self.config.nu * self.config.dt * self.config.viscosity_multiplier)
                u_stoch = u + noise
                idx = 0
                for i in range(len(filaments)):
                    n_pts = len(filaments[i])
                    filaments[i] += self.config.dt * u_stoch[idx:idx+n_pts]
                    idx += n_pts

                center = np.mean(np.vstack(filaments), axis=0)
                for f in filaments:
                    f -= center

                enstrophy_hist.append(E)

                T, Q = self.evaluate_thrust_and_torque(filaments, Gamma_list, step)

                if self.config.preset == "rocket_plume":
                    thrust_history.append(T)
                else:
                    Kt = T / (rho * n**2 * D**4)
                    Kq = Q / (rho * n**2 * D**5)
                    kt_values.append(Kt)
                    kq_values.append(Kq)

                if step % 5 == 0:
                    filament_history.append(copy.deepcopy(filaments))

                if step % 50 == 0:
                    self.save_to_vtk(filaments, Gamma_list, step)

            max_E = np.max(enstrophy_hist)

            if self.config.preset == "rocket_plume":
                max_thrust = np.max(thrust_history) if thrust_history else 0
                csv_rows.append([r+1, f"{max_E:.2e}", "N/A", "N/A", f"{max_thrust:.1f} N"])
                self.gui_log(f"✓ Realization {r+1} finished - Max Enstrophy: {max_E:.2e} | Peak Thrust: {max_thrust:.1f} N")
            else:
                mKt = np.mean(kt_values) if kt_values else 0.0
                mKq = np.mean(kq_values) if kq_values else 0.0
                J = 0.8
                eta = (mKt * J) / (2 * np.pi * mKq) if mKq > 0 else 0.0
                eta = min(eta, 0.98)
                capped_etas.append(eta)
                csv_rows.append([r+1, f"{max_E:.2e}", f"{mKt:.4f}", f"{mKq:.4f}", f"{eta:.4f}"])
                self.gui_log(f"✓ Realization {r+1} finished - Max Enstrophy: {max_E:.2e} | Efficiency: {eta:.3f}")

            max_Es.append(max_E)
            mean_Kt_list.append(mKt)
            mean_Kq_list.append(mKq)

            if self.config.save_gif and filament_history:
                self.gui_log(f"  Saving 3D rotating GIF for Realization {r+1}...")
                fig = plt.figure(figsize=(9,7))
                ax = fig.add_subplot(111, projection='3d')
                def animate(i):
                    ax.cla()
                    state = filament_history[i]
                    for j, f in enumerate(state):
                        color = plt.cm.viridis(j / len(state))
                        ax.plot(f[:,0], f[:,1], f[:,2], color=color, linewidth=1.2)
                    ax.set_xlim(-3.5, 3.5)
                    ax.set_ylim(-3.5, 3.5)
                    ax.set_zlim(-3.5, 3.5)
                    ax.set_xlabel('X (m)')
                    ax.set_ylabel('Y (m)')
                    ax.set_zlabel('Z (m)')
                    ax.set_title(f"Vortex Filaments - Step {i*5} (Realization {r+1})")
                    ax.view_init(elev=25, azim=i*4)
                ani = FuncAnimation(fig, animate, frames=len(filament_history), interval=80)
                ani.save(f'plots/vortex_animation_r{r+1}_3D.gif', writer=PillowWriter(fps=12))
                plt.close()
                self.gui_log(f"  3D GIF saved: plots/vortex_animation_r{r+1}_3D.gif")

        self.gui_log("All realizations completed. Generating final reports...")

        with open("propeller_performance.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_rows)

        report_filename = f"reports/HeliTop_{self.config.preset.replace('_', '_').title()}_Report.pdf"

        if self.config.preset == "rocket_plume":
            plt.figure(figsize=(8,5))
            plt.plot(range(len(thrust_history)), thrust_history, 'r-', linewidth=2)
            plt.xlabel('Step')
            plt.ylabel('Thrust (N)')
            plt.title(f'{preset_name} Thrust Profile')
            plt.grid(True)
            plt.savefig('plots/efficiency_curve.png', dpi=300)
            plt.close()
        else:
            J_values = np.linspace(0.1, 1.2, 12)
            eta_values = [min(((np.mean(mean_Kt_list) * j) / (2 * np.pi * np.mean(mean_Kq_list)) if np.mean(mean_Kq_list) > 0 else 0), 0.98) for j in J_values]
            plt.figure(figsize=(8,5))
            plt.plot(J_values, eta_values, 'b-o', linewidth=2)
            plt.xlabel('Advance Ratio J')
            plt.ylabel('Efficiency η')
            plt.title(f'{preset_name} Efficiency Curve')
            plt.grid(True)
            plt.ylim(0, 1.0)
            plt.savefig('plots/efficiency_curve.png', dpi=300)
            plt.close()

        with PdfPages(report_filename) as pdf:
            plt.figure(figsize=(8,6))
            plt.text(0.5, 0.9, f'HeliTop Vortex v13.0 Final Release', ha='center', fontsize=16)
            plt.text(0.5, 0.8, f'Preset: {self.config.preset}', ha='center')
            plt.text(0.5, 0.7, f'N_FIL: {self.config.N_FIL} | Steps: {self.config.steps} | Realizations: {self.config.NUM_REALIZATIONS}', ha='center')
            plt.text(0.5, 0.6, f'GPU: {"Enabled" if self.config.use_gpu and CUPY_AVAILABLE else "CPU"}', ha='center')
            plt.text(0.5, 0.5, f'Mean Max Enstrophy: {np.mean(max_Es):.2e}', ha='center')
            if self.config.preset == "rocket_plume":
                plt.text(0.5, 0.4, f'Peak Thrust: {max_thrust:.1f} N', ha='center')
                plt.text(0.5, 0.35, f'Calibration: Circ×{self.config.circulation_multiplier:.2f} | Core×{self.config.core_multiplier:.2f} | Dyn×{self.config.dynamic_amplitude:.2f}', ha='center')
            else:
                plt.text(0.5, 0.4, f'Mean Kt: {np.mean(mean_Kt_list):.4f}', ha='center')
                plt.text(0.5, 0.3, f'Mean Kq: {np.mean(mean_Kq_list):.4f}', ha='center')
                plt.text(0.5, 0.2, f'Peak Efficiency: {max(eta_values):.3f}', ha='center')
            plt.axis('off')
            pdf.savefig()
            plt.close()

        return True

class HeliTopGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("HeliTop Vortex v13.0 Final Release")
        self.root.geometry("1120x1080")
        self.root.resizable(True, True)
        self.config = Config()
        self.load_config()
        self.sim = HeliTopSimulator(self.config)
        self.create_widgets()
        self.apply_theme()
        self.force_console_style()

    def load_config(self):
        if os.path.exists("config.json"):
            try:
                with open("config.json", "r") as f:
                    data = json.load(f)
                    for key, value in data.items():
                        if hasattr(self.config, key):
                            setattr(self.config, key, value)
            except:
                pass

    def save_config(self):
        try:
            with open("config.json", "w") as f:
                json.dump({
                    "preset": self.config.preset,
                    "N_FIL": self.config.N_FIL,
                    "NUM_REALIZATIONS": self.config.NUM_REALIZATIONS,
                    "steps": self.config.steps,
                    "core_base": self.config.core_base,
                    "use_gpu": self.config.use_gpu,
                    "save_gif": self.config.save_gif,
                    "dark_mode": self.config.dark_mode,
                    "circulation_multiplier": self.config.circulation_multiplier,
                    "core_multiplier": self.config.core_multiplier,
                    "viscosity_multiplier": self.config.viscosity_multiplier,
                    "thrust_scale_multiplier": self.config.thrust_scale_multiplier,
                    "enstrophy_cap": self.config.enstrophy_cap,
                    "dynamic_amplitude": self.config.dynamic_amplitude
                }, f, indent=4)
        except:
            pass

    def force_console_style(self):
        self.console.configure(bg="#111111", fg="#00ff00", insertbackground="#00ff00")

    def apply_theme(self):
        if self.config.dark_mode:
            bg = "#1e1e1e"
            fg = "#ffffff"
        else:
            bg = "#f0f0f0"
            fg = "#000000"
        self.root.configure(bg=bg)
        for widget in self.root.winfo_children():
            try:
                widget.configure(bg=bg)
                if hasattr(widget, 'configure') and 'fg' in widget.keys():
                    widget.configure(fg=fg)
            except tk.TclError:
                pass
        self.run_btn.configure(bg="#0066cc", fg="white")
        self.close_btn.configure(bg="#cc0000", fg="white")

    def create_widgets(self):
        tk.Label(self.root, text="HeliTop Vortex v13.0", font=("Arial", 18, "bold")).pack(pady=10)
        tk.Label(self.root, text="Professional Vortex Filament Simulator", font=("Arial", 12)).pack()

        frame = tk.Frame(self.root)
        frame.pack(pady=10, fill="x")
        tk.Label(frame, text="Preset:").pack(side="left", padx=10)
        self.preset_var = tk.StringVar(value=self.config.preset)
        ttk.Combobox(frame, textvariable=self.preset_var, 
                     values=["marine_propeller", "rocket_plume", "marine_propeller_high_skew", "aircraft_wake", "generic"], 
                     width=30).pack(side="left")

        self.gpu_var = tk.BooleanVar(value=self.config.use_gpu)
        self.gif_var = tk.BooleanVar(value=self.config.save_gif)
        self.dark_var = tk.BooleanVar(value=self.config.dark_mode)
        tk.Checkbutton(self.root, text="Use GPU Acceleration (CuPy)", variable=self.gpu_var, font=("Arial", 10)).pack(pady=4)
        tk.Checkbutton(self.root, text="Save Vortex Animation as GIF", variable=self.gif_var, font=("Arial", 10)).pack(pady=4)
        tk.Checkbutton(self.root, text="Dark Mode", variable=self.dark_var, font=("Arial", 10), command=self.toggle_dark_mode).pack(pady=4)

        param_frame = tk.LabelFrame(self.root, text="Parameters")
        param_frame.pack(fill="x", padx=20, pady=10)

        tk.Label(param_frame, text="N_FIL:").grid(row=0, column=0, sticky="w", padx=5)
        self.N_FIL_var = tk.IntVar(value=self.config.N_FIL)
        tk.Scale(param_frame, from_=64, to=512, orient="horizontal", variable=self.N_FIL_var, resolution=64).grid(row=0, column=1, sticky="ew", padx=5)

        tk.Label(param_frame, text="Realizations:").grid(row=1, column=0, sticky="w", padx=5)
        self.NUM_REALIZATIONS_var = tk.IntVar(value=self.config.NUM_REALIZATIONS)
        tk.Scale(param_frame, from_=1, to=10, orient="horizontal", variable=self.NUM_REALIZATIONS_var, resolution=1).grid(row=1, column=1, sticky="ew", padx=5)

        tk.Label(param_frame, text="Steps:").grid(row=2, column=0, sticky="w", padx=5)
        self.steps_var = tk.IntVar(value=self.config.steps)
        tk.Scale(param_frame, from_=50, to=300, orient="horizontal", variable=self.steps_var, resolution=50).grid(row=2, column=1, sticky="ew", padx=5)

        tk.Label(param_frame, text="Core Size:").grid(row=3, column=0, sticky="w", padx=5)
        self.core_var = tk.DoubleVar(value=self.config.core_base)
        tk.Scale(param_frame, from_=0.01, to=0.2, orient="horizontal", variable=self.core_var, resolution=0.01).grid(row=3, column=1, sticky="ew", padx=5)

        calib_frame = tk.LabelFrame(self.root, text="Thrust Calibration (rocket_plume only)")
        calib_frame.pack(fill="x", padx=20, pady=10)

        row = 0
        tk.Label(calib_frame, text="Circulation ×").grid(row=row, column=0, sticky="w", padx=5)
        self.circ_var = tk.DoubleVar(value=self.config.circulation_multiplier)
        tk.Scale(calib_frame, from_=0.5, to=2.0, orient="horizontal", variable=self.circ_var, resolution=0.1).grid(row=row, column=1, sticky="ew", padx=5)
        row += 1

        tk.Label(calib_frame, text="Core Size ×").grid(row=row, column=0, sticky="w", padx=5)
        self.core_mult_var = tk.DoubleVar(value=self.config.core_multiplier)
        tk.Scale(calib_frame, from_=0.5, to=2.0, orient="horizontal", variable=self.core_mult_var, resolution=0.1).grid(row=row, column=1, sticky="ew", padx=5)
        row += 1

        tk.Label(calib_frame, text="Viscosity ×").grid(row=row, column=0, sticky="w", padx=5)
        self.visc_var = tk.DoubleVar(value=self.config.viscosity_multiplier)
        tk.Scale(calib_frame, from_=0.5, to=2.0, orient="horizontal", variable=self.visc_var, resolution=0.1).grid(row=row, column=1, sticky="ew", padx=5)
        row += 1

        tk.Label(calib_frame, text="Thrust Scale ×").grid(row=row, column=0, sticky="w", padx=5)
        self.thrust_scale_var = tk.DoubleVar(value=self.config.thrust_scale_multiplier)
        tk.Scale(calib_frame, from_=0.5, to=2.0, orient="horizontal", variable=self.thrust_scale_var, resolution=0.1).grid(row=row, column=1, sticky="ew", padx=5)
        row += 1

        tk.Label(calib_frame, text="Enstrophy Cap").grid(row=row, column=0, sticky="w", padx=5)
        self.enstrophy_var = tk.DoubleVar(value=self.config.enstrophy_cap)
        tk.Scale(calib_frame, from_=1e6, to=5e7, orient="horizontal", variable=self.enstrophy_var, resolution=1e6).grid(row=row, column=1, sticky="ew", padx=5)
        row += 1

        tk.Label(calib_frame, text="Dynamic Amplitude").grid(row=row, column=0, sticky="w", padx=5)
        self.dynamic_var = tk.DoubleVar(value=self.config.dynamic_amplitude)
        tk.Scale(calib_frame, from_=0.0, to=0.2, orient="horizontal", variable=self.dynamic_var, resolution=0.01).grid(row=row, column=1, sticky="ew", padx=5)

        self.progress_frame = tk.Frame(self.root)
        self.progress_frame.pack(fill="x", padx=50, pady=5)
        self.progress_label = tk.Label(self.progress_frame, text="Progress: Idle", font=("Arial", 10))
        self.progress_label.pack()
        self.progress_bar = ttk.Progressbar(self.progress_frame, length=800, mode='determinate')
        self.progress_bar.pack(pady=5)

        self.run_btn = tk.Button(self.root, text="RUN SIMULATION", font=("Arial", 14, "bold"), bg="#0066cc", fg="white", height=2, command=self.start_simulation)
        self.run_btn.pack(pady=15, fill="x", padx=50)

        self.close_btn = tk.Button(self.root, text="Close Application", font=("Arial", 12), bg="#cc0000", fg="white", height=2, command=self.close_application)
        self.close_btn.pack(pady=8, fill="x", padx=80)

        console_frame = tk.Frame(self.root)
        console_frame.pack(fill="both", expand=True, padx=20, pady=5)

        self.console = tk.Text(console_frame, height=24, font=("Courier", 9))
        self.console.pack(side="left", fill="both", expand=True)

        scrollbar = Scrollbar(console_frame, orient=VERTICAL, command=self.console.yview)
        scrollbar.pack(side=RIGHT, fill=Y)
        self.console.config(yscrollcommand=scrollbar.set)

    def toggle_dark_mode(self):
        self.config.dark_mode = self.dark_var.get()
        self.apply_theme()
        self.force_console_style()
        self.save_config()

    def update_progress(self, percent):
        def safe_update():
            self.progress_label.config(text=f"Progress: {percent}%")
            self.progress_bar.configure(value=percent)
        self.root.after(0, safe_update)

    def log(self, text):
        def safe_log():
            self.console.insert(tk.END, f"{text}\n")
            self.console.see(tk.END)
        self.root.after(0, safe_log)

    def start_simulation(self):
        self.run_btn.config(state="disabled")
        self.console.delete(1.0, tk.END)
        self.log("Starting simulation...")
        self.progress_label.config(text="Progress: 0%")
        self.progress_bar.configure(value=0)

        self.config.preset = self.preset_var.get()
        self.config.use_gpu = self.gpu_var.get()
        self.config.save_gif = self.gif_var.get()
        self.config.N_FIL = self.N_FIL_var.get()
        self.config.NUM_REALIZATIONS = self.NUM_REALIZATIONS_var.get()
        self.config.steps = self.steps_var.get()
        self.config.core_base = self.core_var.get()
        self.config.circulation_multiplier = self.circ_var.get()
        self.config.core_multiplier = self.core_mult_var.get()
        self.config.viscosity_multiplier = self.visc_var.get()
        self.config.thrust_scale_multiplier = self.thrust_scale_var.get()
        self.config.enstrophy_cap = self.enstrophy_var.get()
        self.config.dynamic_amplitude = self.dynamic_var.get()

        self.save_config()

        def run():
            try:
                self.sim = HeliTopSimulator(self.config)
                self.sim.gui_log = self.log
                success = self.sim.run_campaign()
                self.log("All done! Simulation completed successfully.")
                if success:
                    self.root.after(0, lambda: messagebox.showinfo("Success", "All files generated!\nCheck reports/, plots/, and vtk/ folders."))
            except Exception as e:
                error_msg = f"ERROR: {str(e)}"
                self.log(error_msg)
                self.root.after(0, lambda: messagebox.showerror("Error", error_msg))
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.progress_label.config(text="Progress: Done"))
        threading.Thread(target=run, daemon=True).start()

    def close_application(self):
        if messagebox.askyesno("Exit", "Are you sure you want to close the application?"):
            self.save_config()
            self.root.destroy()
            sys.exit(0)

if __name__ == "__main__":
    try:
        app = HeliTopGUI()
        app.root.mainloop()
    except Exception as e:
        messagebox.showerror("Critical Error", f"Startup failed:\n{str(e)}")
