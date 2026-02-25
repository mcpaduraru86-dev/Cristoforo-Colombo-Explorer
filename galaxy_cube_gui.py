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

# 2. THE DEEP LOADER (STRIPS 'Z' ERRORS)
print("--- NORTH STAR V7.1: LOCKING TARGET ---")
try:
    with fits.open(FULL_PATH) as hdul:
        data = hdul[1].data
        raw_cols = data.names
        # Create a dictionary: {CLEAN_UPPER_NAME: ORIGINAL_NAME}
        col_map = {str(c).strip().upper(): c for c in raw_cols}
        
        # Priority search for the Redshift column
        z_col = None
        for candidate in ['Z', 'Z_HP', 'REDSHIFT', 'Z_NOT4CLUS', 'Z_TARGET', 'Z_FLUX']:
            if candidate in col_map:
                z_col = col_map[candidate]
                break
        
        # If still not found, grab the first column that has a 'Z' in it
        if not z_col:
            z_col = next((col_map[c] for c in col_map if 'Z' in c), None)

        if not z_col:
            raise KeyError(f"Redshift column not found. Available: {list(col_map.keys())}")

        ra_col = col_map.get('RA', 'RA')
        dec_col = col_map.get('DEC', 'DEC')

        mask = np.isfinite(data[z_col]) & (data[z_col] > 0.001)
        dist = cosmo.comoving_distance(data[z_col][mask]).value
        ra_r, de_r = np.deg2rad(data[ra_col][mask]), np.deg2rad(data[dec_col][mask])
        coords_3d = np.vstack([dist*np.cos(de_r)*np.cos(ra_r), dist*np.cos(de_r)*np.sin(ra_r), dist*np.sin(de_r)]).T
        tree = cKDTree(coords_3d)
        print(f"Engine Ready. Using column: {z_col}")
except Exception as e:
    print(f"FATAL LOADER ERROR: {e}")
    os._exit(1)

# 3. NORTH STAR CLIMBING LOGIC
def build_squeezed_cube(center, axis, size, squeeze):
    u = normalize(axis)
    temp = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    v = normalize(np.cross(u, temp)); w = np.cross(u, v)
    return np.array([center + a*(size*u) + b*(size*squeeze*v) + c*(size*squeeze*w) for a in [0,1] for b in [0,1] for c in [0,1]])

def run_north_star_search():
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, ">>> TRACKING THE NORTH STAR (CLIMBING TO 80%+)...\n")
    
    r_f, d_f, z_f = float(entry_ra.get()), float(entry_dec.get()), float(entry_z.get())
    dist_t = cosmo.comoving_distance(z_f).value
    ra_r, de_r = np.deg2rad(r_f), np.deg2rad(d_f)
    current_pt = np.array([dist_t*np.cos(de_r)*np.cos(ra_r), dist_t*np.cos(de_r)*np.sin(ra_r), dist_t*np.sin(de_r)])
    current_best_rate = 0.0

    for i in range(30): # Climbing steps
        found_step = False
        # Search a small cloud around the current point
        search_cloud = (np.random.rand(50, 3) - 0.5) * 10.0 
        
        for offset in search_cloud:
            test_pt = current_pt + offset
            vtxs = build_squeezed_cube(test_pt, BASE_AXIS, EDGE_MPC, BEST_SQUEEZE)
            hits = sum(1 for v in vtxs if tree.query(v)[0] < DENSITY_RADIUS_MPC)
            rate = (hits / 8.0) * 100
            
            if rate > current_best_rate:
                current_best_rate = rate
                current_pt = test_pt
                found_step = True
                text_output.insert(tk.END, f"Climb {i+1}: {current_best_rate:.1f}% hit found!\n")
                root.update_idletasks()
                if current_best_rate >= 87.5: break
        
        if not found_step: break

    # Final Coordinates
    final_dist = np.linalg.norm(current_pt)
    final_dec = np.rad2deg(np.arcsin(current_pt[2] / final_dist))
    final_ra = np.rad2deg(np.arctan2(current_pt[1], current_pt[0]))
    if final_ra < 0: final_ra += 360
    
    text_output.insert(tk.END, f"\n🏆 CAPITAL NODE LOCATED 🏆\nRA: {final_ra:.4f}\nDEC: {final_dec:.4f}\nRESONANCE: {current_best_rate:.1f}%")

# 4. DASHBOARD
root = tk.Tk(); root.title("V7.1 - North Star Tracker")
entry_frame = tk.Frame(root); entry_frame.pack(pady=10)
entry_ra = tk.Entry(entry_frame, width=10); entry_ra.insert(0, "189.481"); entry_ra.pack(side="left")
entry_dec = tk.Entry(entry_frame, width=10); entry_dec.insert(0, "31.969"); entry_dec.pack(side="left")
entry_z = tk.Entry(entry_frame, width=10); entry_z.insert(0, "0.024"); entry_z.pack(side="left")
text_output = tk.Text(root, width=75, height=20); text_output.pack(padx=10, pady=10)
tk.Button(root, text="🌟 TRACK NORTH STAR", bg="cyan", font=('bold'), command=lambda: threading.Thread(target=run_north_star_search, daemon=True).start()).pack(pady=10)
root.mainloop()