import os
import threading
import numpy as np
import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
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

FULL_PATH = r"C:\Users\Marian Catalin\Documents\Ianus Paduraru\BGS_BRIGHT_full.dat.fits"

# 2. THE ENGINE (The "Detective" Loader)
print("Engaging Indexers... Searching for data columns.")
with fits.open(FULL_PATH) as hdul:
    data = hdul[1].data
    cols = data.names
    print(f"I found these columns in your file: {cols}")
    
    # Logic to find Redshift even if it has a strange name
    z_col = None
    for candidate in ['Z', 'z', 'Z_HP', 'Z_not4clus', 'TARGET_Z', 'REDSHIFT']:
        if candidate in cols:
            z_col = candidate
            break
    
    if not z_col:
        # Emergency backup: just take the first column that starts with Z
        z_col = next((c for c in cols if c.upper().startswith('Z')), None)

    if not z_col:
        print("CRITICAL ERROR: I can't find a Redshift column. Please check the list above.")
        os._exit(1)

    print(f"Using '{z_col}' as our Redshift tracker.")
    mask = np.isfinite(data[z_col]) & (data[z_col] > 0.001)
    ra, dec, z = data["RA"][mask], data["DEC"][mask], data[z_col][mask]
    targetid = data["TARGETID"][mask].astype(np.int64)

# 3D Mapping
Dc = cosmo.comoving_distance(z).value
ra_rad, dec_rad = np.deg2rad(ra), np.deg2rad(dec)
coords = np.vstack([
    Dc * np.cos(dec_rad) * np.cos(ra_rad),
    Dc * np.cos(dec_rad) * np.sin(ra_rad),
    Dc * np.sin(dec_rad)
]).T

tree = cKDTree(coords) 
id_to_index = {int(tid): i for i, tid in enumerate(targetid)}
last_corner_data = None
last_target_id = None

# 3. THE TILT LOGIC
def build_cube(center, axis, size):
    u = normalize(axis)
    temp = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    v = normalize(np.cross(u, temp))
    w = np.cross(u, v)
    e1, e2, e3 = size * u, size * v, size * w
    return np.array([center + a*e1 + b*e2 + c*e3 for a in [0,1] for b in [0,1] for c in [0,1]])

def run_tilt_scan():
    global BASE_AXIS
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, f">>> PORT ROYALE PRECISION TILT\nLocked Scale: {EDGE_MPC} Mpc\n")
    
    best_rate = 0
    best_vector = BASE_AXIS
    
    for i in range(100):
        nudge = (np.random.rand(3) - 0.5) * 0.04 # Fine adjustment nudge
        test_vec = normalize(BASE_AXIS + nudge)
        
        sample_indices = np.random.choice(len(coords), 100, replace=False)
        hits = 0
        for idx in sample_indices:
            vtxs = build_cube(coords[idx], test_vec, EDGE_MPC)
            for v in vtxs:
                if tree.query(v)[0] < DENSITY_RADIUS_MPC: hits += 1
        
        rate = (hits / 800.0) * 100
        if rate > best_rate:
            best_rate = rate
            best_vector = test_vec
            text_output.insert(tk.END, f"Signal Improving: {rate:.1f}% | Vector locked.\n")
            root.update_idletasks()
            
    BASE_AXIS = best_vector
    text_output.insert(tk.END, f"\n--- SCAN FINISHED ---\nPeak Port Royale Resonance: {best_rate:.1f}%\n")
    messagebox.showinfo("Target Reached", f"New Peak: {best_rate:.1f}%")

# 4. DASHBOARD
root = tk.Tk(); root.title("Port Royale V2.1 - Column Blind Navigator")

frame_gps = tk.Frame(root); frame_gps.pack(pady=10)
entry_ra = tk.Entry(frame_gps, width=8); entry_ra.insert(0, "189.481"); entry_ra.pack(side="left")
entry_dec = tk.Entry(frame_gps, width=8); entry_dec.insert(0, "31.969"); entry_dec.pack(side="left")
entry_z = tk.Entry(frame_gps, width=8); entry_z.insert(0, "0.024"); entry_z.pack(side="left")

def lock_target():
    global last_corner_data, last_target_id
    try:
        r_f, d_f, z_f = float(entry_ra.get()), float(entry_dec.get()), float(entry_z.get())
        dist_t = cosmo.comoving_distance(z_f).value
        ra_r, de_r = np.deg2rad(r_f), np.deg2rad(d_f)
        target_xyz = np.array([dist_t*np.cos(de_r)*np.cos(ra_r), dist_t*np.cos(de_r)*np.sin(ra_r), dist_t*np.sin(de_r)])
        _, idx = tree.query(target_xyz)
        last_target_id = int(targetid[idx])
        vtxs = build_cube(coords[idx], BASE_AXIS, EDGE_MPC)
        last_corner_data = [{"tid": int(targetid[tree.query(v)[1]]), "lab": "GALAXY" if tree.query(v)[0] < DENSITY_RADIUS_MPC else "VOID"} for v in vtxs]
        text_output.insert(tk.END, f"Target {last_target_id} Locked.\n")
    except: messagebox.showerror("Error", "Check Inputs")

tk.Button(root, text="🎯 LOCK TARGET", command=lock_target, bg="#c1f0c1").pack(pady=5)
text_output = tk.Text(root, width=70, height=15); text_output.pack(padx=10, pady=10)

tk.Button(root, text="🌀 HIGH-PRECISION TILT SCAN", bg="orange", font=('bold'), 
          command=lambda: threading.Thread(target=run_tilt_scan, daemon=True).start()).pack(pady=5)

def show_3d():
    if not last_target_id: return
    fig = plt.figure(figsize=(10,8)); ax = fig.add_subplot(111, projection='3d')
    t_pos = coords[id_to_index[last_target_id]]
    idx_neigh = tree.query_ball_point(t_pos, r=100); neigh = coords[idx_neigh]
    ax.scatter(neigh[:,0], neigh[:,1], neigh[:,2], s=1, alpha=0.1, color='gray')
    ax.scatter(t_pos[0], t_pos[1], t_pos[2], s=120, color='red')
    for c in last_corner_data:
        p = coords[id_to_index[c['tid']]]
        ax.scatter(p[0], p[1], p[2], s=70, color='blue' if c['lab']=='GALAXY' else 'gold')
    plt.show()

tk.Button(root, text="🌌 SHOW 3D STRUCTURE", bg="#3498db", fg="white", command=show_3d).pack(pady=5)
root.mainloop()