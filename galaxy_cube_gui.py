import os
import threading
import numpy as np
import tkinter as tk
from tkinter import messagebox
from astropy.io import fits
from scipy.spatial import cKDTree
from astropy.cosmology import FlatLambdaCDM

# 1. THE STEERING
def normalize(v): 
    norm = np.linalg.norm(v)
    return v / norm if norm > 0 else v

BASE_AXIS = normalize(np.array([0.311, 0.814, 0.490]))
cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089)
DENSITY_RADIUS_MPC = 10.0
EDGE_MPC = 137.5 
BEST_SQUEEZE = 0.950 

FULL_PATH = r"C:\Users\Marian Catalin\Documents\Ianus Paduraru\BGS_BRIGHT_full.dat.fits"

# 2. THE DEEP LOADER (STRIP & FIND)
print("--- MIRROR MASTER V6.1: DEEP SCANNING COLUMNS ---")
try:
    with fits.open(FULL_PATH) as hdul:
        data = hdul[1].data
        # Clean the column names: strip spaces and make uppercase
        raw_cols = data.names
        clean_cols = [str(c).strip().upper() for c in raw_cols]
        
        # Mapping clean names back to original names
        col_map = dict(zip(clean_cols, raw_cols))

        z_key = None
        # Priority list for Redshift
        for candidate in ['Z', 'Z_HP', 'REDSHIFT', 'Z_NOT4CLUS', 'Z_TARGET']:
            if candidate in clean_cols:
                z_key = col_map[candidate]
                break
        
        if not z_key: # Last resort: find any column that is just 'Z' with extra bits
            z_key = next((col_map[c] for c in clean_cols if 'Z' in c), None)

        if not z_key:
            raise KeyError(f"Could not find Redshift. Available: {clean_cols}")

        ra_key = col_map.get('RA', 'RA')
        dec_key = col_map.get('DEC', 'DEC')

        print(f"Locked on: Z='{z_key}', RA='{ra_key}', DEC='{dec_key}'")

        mask = np.isfinite(data[z_key]) & (data[z_key] > 0.001)
        dist = cosmo.comoving_distance(data[z_key][mask]).value
        ra_r, de_r = np.deg2rad(data[ra_key][mask]), np.deg2rad(data[dec_key][mask])
        coords_3d = np.vstack([dist*np.cos(de_r)*np.cos(ra_r), dist*np.cos(de_r)*np.sin(ra_r), dist*np.sin(de_r)]).T
        tree = cKDTree(coords_3d)
except Exception as e:
    print(f"FATAL ERROR: {e}")
    os._exit(1)

# 3. PIVOT LOGIC
def build_squeezed_cube(center, axis, size, squeeze):
    u = normalize(axis)
    temp = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    v = normalize(np.cross(u, temp)); w = np.cross(u, v)
    return np.array([center + a*(size*u) + b*(size*squeeze*v) + c*(size*squeeze*w) for a in [0,1] for b in [0,1] for c in [0,1]])

def run_pivot_scan():
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, ">>> MIRROR PIVOT: SEEKING THE HIDDEN CURVE...\n")
    
    try:
        r_f, d_f, z_f = float(entry_ra.get()), float(entry_dec.get()), float(entry_z.get())
        dist_t = cosmo.comoving_distance(z_f).value
        ra_r, de_r = np.deg2rad(r_f), np.deg2rad(d_f)
        center_pt = np.array([dist_t*np.cos(de_r)*np.cos(ra_r), dist_t*np.cos(de_r)*np.sin(ra_r), dist_t*np.sin(de_r)])

        for step in range(-1, -4, -1):
            jump = center_pt + (step * EDGE_MPC * BASE_AXIS)
            best_rate = 0
            
            for _ in range(120): # Higher precision wiggle
                nudge = (np.random.rand(3) - 0.5) * 0.20
                test_axis = normalize(BASE_AXIS + nudge)
                _, neighbors = tree.query(jump, k=100)
                hits = sum(1 for idx in neighbors for v in build_squeezed_cube(coords_3d[idx], test_axis, EDGE_MPC, BEST_SQUEEZE) if tree.query(v)[0] < DENSITY_RADIUS_MPC)
                rate = (hits / 800.0) * 100
                if rate > best_rate: best_rate = rate
            
            text_output.insert(tk.END, f"Mirror Step {abs(step)} (-{abs(int(step*EDGE_MPC))} Mpc): Peak Found at {best_rate:.1f}%\n")
            root.update_idletasks()
        text_output.insert(tk.END, "\n--- IMPERIAL CURVATURE MAPPED ---\n")
    except Exception as e:
        text_output.insert(tk.END, f"\nError: {e}")

# 4. DASHBOARD
root = tk.Tk(); root.title("V6.1 - Pivot Deep Scan")
entry_frame = tk.Frame(root); entry_frame.pack(pady=10)
entry_ra = tk.Entry(entry_frame, width=10); entry_ra.insert(0, "189.481"); entry_ra.pack(side="left")
entry_dec = tk.Entry(entry_frame, width=10); entry_dec.insert(0, "31.969"); entry_dec.pack(side="left")
entry_z = tk.Entry(entry_frame, width=10); entry_z.insert(0, "0.024"); entry_z.pack(side="left")

text_output = tk.Text(root, width=75, height=20); text_output.pack(padx=10, pady=10)
tk.Button(root, text="🌀 FIRE PIVOT SCAN", bg="orange", font=('bold'), 
          command=lambda: threading.Thread(target=run_pivot_scan, daemon=True).start()).pack(pady=10)
root.mainloop()