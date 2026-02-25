import os
import threading
import numpy as np
import tkinter as tk
from tkinter import messagebox
from astropy.io import fits
from scipy.spatial import cKDTree
from astropy.cosmology import FlatLambdaCDM

# 1. THE STEERING (LOCKED WINNER SPECS)
def normalize(v): 
    norm = np.linalg.norm(v)
    return v / norm if norm > 0 else v

BASE_AXIS = normalize(np.array([0.311, 0.814, 0.490]))
cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089)
DENSITY_RADIUS_MPC = 10.0
EDGE_MPC = 137.5 
BEST_SQUEEZE = 0.950 

FULL_PATH = r"C:\Users\Marian Catalin\Documents\Ianus Paduraru\BGS_BRIGHT_full.dat.fits"

# 2. THE BULLETPROOF LOADER
print("--- MIRROR MASTER: CHECKING DATA MANIFEST ---")
try:
    with fits.open(FULL_PATH) as hdul:
        data = hdul[1].data
        cols = data.names
        print(f"I see these columns: {cols}")
        
        # Detect Redshift (Z)
        z_col = next((c for c in cols if c.upper() in ['Z', 'Z_HP', 'REDSHIFT', 'Z_TARGET']), None)
        if not z_col: z_col = next((c for c in cols if 'Z' in c.upper()), None)
        
        # Detect RA and DEC
        ra_col = next((c for c in cols if c.upper() in ['RA', 'TARGET_RA']), 'RA')
        dec_col = next((c for c in cols if c.upper() in ['DEC', 'TARGET_DEC']), 'DEC')
        
        if not z_col: raise KeyError("No Redshift column found.")
        print(f"Locked on: RA={ra_col}, DEC={dec_col}, Z={z_col}")

        mask = np.isfinite(data[z_col]) & (data[z_col] > 0.001)
        z_vals = data[z_col][mask]
        ra_vals = data[ra_col][mask]
        dec_vals = data[dec_col][mask]
        
        dist = cosmo.comoving_distance(z_vals).value
        ra_r, de_r = np.deg2rad(ra_vals), np.deg2rad(dec_vals)
        coords_3d = np.vstack([dist*np.cos(de_r)*np.cos(ra_r), dist*np.cos(de_r)*np.sin(ra_r), dist*np.sin(de_r)]).T
        tree = cKDTree(coords_3d)
except Exception as e:
    print(f"LOADER ERROR: {e}")
    os._exit(1)

# 3. MIRROR LOGIC
def build_squeezed_cube(center, axis, size, squeeze):
    u = normalize(axis)
    temp = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    v = normalize(np.cross(u, temp)); w = np.cross(u, v)
    return np.array([center + a*(size*u) + b*(size*squeeze*v) + c*(size*squeeze*w) for a in [0,1] for b in [0,1] for c in [0,1]])

def run_mirror_scan():
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, ">>> JANUS MIRROR: CHECKING THE REARVIEW...\n")
    
    try:
        r_f, d_f, z_f = float(entry_ra.get()), float(entry_dec.get()), float(entry_z.get())
        dist_t = cosmo.comoving_distance(z_f).value
        ra_r, de_r = np.deg2rad(r_f), np.deg2rad(d_f)
        center_pt = np.array([dist_t*np.cos(de_r)*np.cos(ra_r), dist_t*np.cos(de_r)*np.sin(ra_r), dist_t*np.sin(de_r)])

        for step in range(-1, -4, -1):
            jump = center_pt + (step * EDGE_MPC * BASE_AXIS)
            _, neighbors = tree.query(jump, k=150)
            hits = sum(1 for idx in neighbors for v in build_squeezed_cube(coords_3d[idx], BASE_AXIS, EDGE_MPC, BEST_SQUEEZE) if tree.query(v)[0] < DENSITY_RADIUS_MPC)
            
            rate = (hits / (150 * 8)) * 100
            text_output.insert(tk.END, f"Mirror Step {abs(step)} (-{abs(int(step*EDGE_MPC))} Mpc): {rate:.1f}%\n")
            root.update_idletasks()
        text_output.insert(tk.END, "\n--- SCAN COMPLETE ---\n")
    except Exception as e:
        text_output.insert(tk.END, f"\nError: {e}")

# 4. DASHBOARD
root = tk.Tk(); root.title("Port Royale V5.2 - Mirror Master")
entry_frame = tk.Frame(root); entry_frame.pack(pady=10)
entry_ra = tk.Entry(entry_frame, width=10); entry_ra.insert(0, "189.481"); entry_ra.pack(side="left")
entry_dec = tk.Entry(entry_frame, width=10); entry_dec.insert(0, "31.969"); entry_dec.pack(side="left")
entry_z = tk.Entry(entry_frame, width=10); entry_z.insert(0, "0.024"); entry_z.pack(side="left")

text_output = tk.Text(root, width=75, height=18); text_output.pack(padx=10, pady=10)
tk.Button(root, text="🔄 RUN MIRROR SCAN (BACKWARD)", bg="blue", fg="white", font=('bold'), 
          command=lambda: threading.Thread(target=run_mirror_scan, daemon=True).start()).pack(pady=10)
root.mainloop()