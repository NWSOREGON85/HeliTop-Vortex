import tkinter as tk
from tkinter import ttk, messagebox
import threading
import subprocess
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import csv
from dataclasses import dataclass

# ====================== CORE SIMULATION ENGINE ======================
@dataclass
class Config:
    preset: str = "marine_propeller"
    N_FIL: int = 256
    NUM_REALIZATIONS: int = 10
    steps: int = 200
    dt: float = 0.0015
    core_base: float = 0.08
    nu: float = 0.0005
    d0_recon: float = 0.12
    alpha_recon: float = 0.75
    angle_threshold_deg: float = 35.0

class HeliTopSimulator:
    def __init__(self, config):
        self.config = config

    def get_dl(self, r):
        return np.roll(r, -1, axis=0) - r

    def biot_savart_induced(self, filaments, Gamma_list, core=0.08):
        all_pts = np.vstack(filaments)
        u = np.zeros_like(all_pts, dtype=np.float64)
        for fil, G in zip(filaments, Gamma_list):
            dl = self.get_dl(fil)
            R = all_pts[:, np.newaxis, :] - fil
            R2 = np.sum(R**2, axis=-1)
            R2_safe = np.maximum(R2, core**2)
            cross = np.cross(dl[np.newaxis, :, :], R)
            factor = (G / (4 * np.pi)) * np.sqrt(R2_safe) / (R2_safe + core**2)
            u += np.sum(factor[..., np.newaxis] * cross, axis=1)
        return u

    def enstrophy_proxy(self, filaments, Gamma_list):
        E = 0.0
        for r, G in zip(filaments, Gamma_list):
            E += G**2 * np.sum(np.linalg.norm(self.get_dl(r), axis=1))
        return E

    def adaptive_regrid(self, r, stretch_threshold=1.5):
        if len(r) < 2: return r.copy()
        dr = np.diff(r, axis=0)
        lengths = np.linalg.norm(dr, axis=1)
        if np.max(lengths) <= stretch_threshold: return r.copy()
        new_r = [r[0]]
        for i in range(len(lengths)):
            new_r.append(r[i+1])
            if lengths[i] > stretch_threshold:
                new_r.append(0.5 * (r[i] + r[i+1]))
        return np.array(new_r)

    def reconnect_filaments(self, filaments, Gamma_list, H_top):
        d_recon = self.config.d0_recon / (1 + self.config.alpha_recon * max(H_top, 0))
        new_filaments = [f.copy() for f in filaments]
        new_Gamma = Gamma_list.copy()
        recon_count = 0
        for i in range(len(filaments)):
            for j in range(i+1, len(filaments)):
                r1, r2 = filaments[i], filaments[j]
                dl1, dl2 = self.get_dl(r1), self.get_dl(r2)
                for k in range(len(r1)-1):
                    for m in range(len(r2)-1):
                        dist = np.linalg.norm((r1[k]+r1[k+1])/2 - (r2[m]+r2[m+1])/2)
                        if dist < d_recon:
                            cos_angle = np.dot(dl1[k], dl2[m]) / (np.linalg.norm(dl1[k])*np.linalg.norm(dl2[m])+1e-12)
                            angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
                            if angle < self.config.angle_threshold_deg or angle > 180 - self.config.angle_threshold_deg:
                                tail_i = new_filaments[i][-1].copy()
                                tail_j = new_filaments[j][-1].copy()
                                new_filaments[i][-1] = tail_j
                                new_filaments[j][-1] = tail_i
                                recon_count += 1
                                break
                    else:
                        continue
                    break
        return new_filaments, new_Gamma, recon_count

    def get_preset_data(self):
        if self.config.preset == "marine_propeller":
            filaments = []
            Gamma_list = []
            num_blades = 4
            radius = 2.5
            pitch_ratio = 0.8
            for b in range(num_blades):
                theta = np.linspace(0, 3*np.pi, self.config.N_FIL)
                r = radius * (0.6 + 0.4 * np.sin(theta * 2))
                x = r * np.cos(theta + b * 2 * np.pi / num_blades)
                y = r * np.sin(theta + b * 2 * np.pi / num_blades)
                z = (pitch_ratio * radius * theta / (2 * np.pi)) + 0.1 * np.sin(theta * 3)
                blade = np.stack((x, y, z), axis=1).astype(np.float32)
                filaments.append(blade)
                Gamma_list.append(1.0)
            for b in range(num_blades):
                theta_tip = np.linspace(0, 6*np.pi, self.config.N_FIL // 2)
                r_tip = radius
                x_tip = r_tip * np.cos(theta_tip + b * 2 * np.pi / num_blades + 0.2)
                y_tip = r_tip * np.sin(theta_tip + b * 2 * np.pi / num_blades + 0.2)
                z_tip = pitch_ratio * radius * theta_tip / (2 * np.pi) + 0.3
                tip_vortex = np.stack((x_tip, y_tip, z_tip), axis=1).astype(np.float32)
                filaments.append(tip_vortex)
                Gamma_list.append(0.6)
            hub_theta = np.linspace(0, 8*np.pi, self.config.N_FIL // 3)
            hub = np.stack((0.3 * np.cos(hub_theta), 0.3 * np.sin(hub_theta), 0.5 * hub_theta / np.pi), axis=1)
            filaments.append(hub.astype(np.float32))
            Gamma_list.append(-0.4)
            return filaments, Gamma_list, radius * 2
        else:
            # generic fallback
            filaments = []
            Gamma_list = []
            for i in range(6):
                radius = 0.8 + 0.25 * i
                theta = np.linspace(0, 2*np.pi, self.config.N_FIL)
                z = -0.5 * i
                r = np.stack((radius * np.cos(theta), radius * np.sin(theta), z), axis=1)
                filaments.append(r.astype(np.float32))
                Gamma_list.append(1.0 if i % 2 == 0 else -1.0)
            return filaments, Gamma_list, 5.0

    def evaluate_thrust_and_torque(self, filaments, Gamma_list):
        thrust = 0.0
        torque = 0.0
        for fil, G in zip(filaments, Gamma_list):
            dl = self.get_dl(fil)
            thrust += G * np.sum(dl[:,2])
            torque += G * np.sum(fil[:,0] * dl[:,1] - fil[:,1] * dl[:,0])
        return abs(thrust), abs(torque)

    def save_to_vtk(self, filaments, Gamma_list, step):
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
        max_Es = []
        mean_Kt = []
        mean_Kq = []
        csv_rows = [["Realization", "Max_Enstrophy", "Mean_Kt", "Mean_Kq", "Mean_Efficiency"]]
        rho = 1025.0
        n = 10.0
        filaments, Gamma_list, diameter = self.get_preset_data()
        D = diameter
        for r in range(self.config.NUM_REALIZATIONS):
            enstrophy_hist = []
            kt_values = []
            kq_values = []
            for step in range(self.config.steps):
                filaments = [self.adaptive_regrid(f) for f in filaments]
                E = self.enstrophy_proxy(filaments, Gamma_list)
                filaments, Gamma_list, _ = self.reconnect_filaments(filaments, Gamma_list, 0.0)
                u = self.biot_savart_induced(filaments, Gamma_list, self.config.core_base)
                noise = np.random.randn(*u.shape) * np.sqrt(2 * self.config.nu * 1.0 * self.config.dt)
                u_stoch = u + noise
                idx = 0
                for i in range(len(filaments)):
                    n_pts = len(filaments[i])
                    filaments[i] += self.config.dt * u_stoch[idx:idx+n_pts]
                    idx += n_pts
                enstrophy_hist.append(E)
                if step % 50 == 0:
                    T, Q = self.evaluate_thrust_and_torque(filaments, Gamma_list)
                    Kt = T / (rho * n**2 * D**4)
                    Kq = Q / (rho * n**2 * D**5)
                    kt_values.append(Kt)
                    kq_values.append(Kq)
                    self.save_to_vtk(filaments, Gamma_list, step)
            max_E = np.max(enstrophy_hist)
            mKt = np.mean(kt_values)
            mKq = np.mean(kq_values)
            J = 0.8
            eta = (mKt * J) / (2 * np.pi * mKq) if mKq > 0 else 0
            max_Es.append(max_E)
            mean_Kt.append(mKt)
            mean_Kq.append(mKq)
            csv_rows.append([r+1, f"{max_E:.4f}", f"{mKt:.4f}", f"{mKq:.4f}", f"{eta:.4f}"])
        # Save CSV
        with open("propeller_performance.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(csv_rows)
        # Efficiency curve
        J_values = np.linspace(0.1, 1.2, 12)
        eta_values = [(np.mean(mean_Kt) * j) / (2 * np.pi * np.mean(mean_Kq)) for j in J_values]
        plt.figure(figsize=(8,5))
        plt.plot(J_values, eta_values, 'b-o', linewidth=2)
        plt.xlabel('Advance Ratio J')
        plt.ylabel('Efficiency η')
        plt.title('Propeller Efficiency Curve')
        plt.grid(True)
        plt.savefig('plots/efficiency_curve.png', dpi=300)
        plt.close()
        # PDF Report
        with PdfPages('reports/HeliTop_Propeller_Report.pdf') as pdf:
            plt.figure(figsize=(8,6))
            plt.text(0.5, 0.9, 'HeliTop Vortex v8.0 Report', ha='center', fontsize=16)
            plt.text(0.5, 0.7, f'Preset: {self.config.preset}', ha='center')
            plt.text(0.5, 0.5, f'Mean Max Enstrophy: {np.mean(max_Es):.2f}', ha='center')
            plt.text(0.5, 0.4, f'Mean Kt: {np.mean(mean_Kt):.4f}', ha='center')
            plt.text(0.5, 0.3, f'Mean Kq: {np.mean(mean_Kq):.4f}', ha='center')
            plt.text(0.5, 0.2, f'Peak Efficiency: {max(eta_values):.3f}', ha='center')
            plt.axis('off')
            pdf.savefig()
            plt.close()
        return True

# ====================== GUI ======================
class HeliTopGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("HeliTop Vortex v8.0")
        self.root.geometry("900x700")
        self.config = Config()
        self.sim = HeliTopSimulator(self.config)
        self.create_widgets()

    def create_widgets(self):
        tk.Label(self.root, text="HeliTop Vortex v8.0", font=("Arial", 18, "bold")).pack(pady=10)
        tk.Label(self.root, text="Professional Vortex Filament Simulator", font=("Arial", 12)).pack()

        frame = tk.Frame(self.root)
        frame.pack(pady=10, fill="x")
        tk.Label(frame, text="Preset:").pack(side="left", padx=10)
        self.preset_var = tk.StringVar(value=self.config.preset)
        ttk.Combobox(frame, textvariable=self.preset_var, values=["marine_propeller", "aircraft_wake", "generic"], width=20).pack(side="left")

        self.run_btn = tk.Button(self.root, text="RUN SIMULATION", font=("Arial", 14, "bold"), bg="#0066cc", fg="white", height=2, command=self.start_simulation)
        self.run_btn.pack(pady=20, fill="x", padx=50)

        self.console = tk.Text(self.root, height=15, bg="#111111", fg="#00ff00", font=("Courier", 9))
        self.console.pack(fill="both", padx=20, pady=5)

    def log(self, text):
        self.console.insert(tk.END, text + "\n")
        self.console.see(tk.END)
        self.root.update()

    def start_simulation(self):
        self.run_btn.config(state="disabled")
        self.console.delete(1.0, tk.END)
        self.log("Starting simulation...")

        def run():
            self.config.preset = self.preset_var.get()
            self.sim = HeliTopSimulator(self.config)
            success = self.sim.run_campaign()
            self.log("Simulation complete!")
            if success:
                messagebox.showinfo("Success", "All files generated!\nCheck reports/, plots/, and vtk/ folders.")
            self.run_btn.config(state="normal")

        threading.Thread(target=run, daemon=True).start()

if __name__ == "__main__":
    HeliTopGUI()
