import os
import threading
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
from astropy.io import fits
from scipy.spatial import cKDTree
from astropy.cosmology import FlatLambdaCDM

# ============================================================
# 1. SETUP & COSMOLOGY
# ============================================================
def build_cosmology():
    return FlatLambdaCDM(H0=67.74, Om0=0.3089)

def radec_z_to_xyz(ra, dec, z, cosmo):
    Dc = cosmo.comoving_distance(z).value
    ra_rad, dec_rad = np.deg2rad(ra), np.deg2rad(dec)
    x = Dc * np.cos(dec_rad) * np.cos(ra_rad)
    y = Dc * np.cos(dec_rad) * np.sin(ra_rad)
    z3d = Dc * np.sin(dec_rad)
    return np.array([x, y, z3d])

# Paths - Double check these match your computer!
FULL_PATH = r"C:\Users\Marian Catalin\Documents\Ianus Paduraru\BGS_BRIGHT_full.dat.fits"
CLUST_PATH = r"C:\Users\Marian Catalin\Documents\Ianus Paduraru\BGS_BRIGHT_NGC_clustering.dat.fits"

DENSITY_RADIUS_MPC = 10.0
EDGE_MPC = 137.5
BASE_AXIS = np.array([0.303, 0.808, 0.505], dtype=float)
BASE_AXIS /= np.linalg.norm(BASE_AXIS)

# ============================================================
# 2. DATA LOADING (The 4.2 Million Haul)
# ============================================================
print("Loading 5GB Catalog... Please wait.")
with fits.open(FULL_PATH) as hdul:
    full_data = hdul[1].data
    full_tid = full_data["TARGETID"].astype(np.int64)
    # Check for correct Redshift column
    possible_z = ['Z', 'Z_not4clus', 'Z_HP', 'Z_PHOT']
    z_col = next((n for n in possible_z if n in full_data.names), 'Z')

# Load the clustering Z map for accuracy
clust_z_map = {}
if os.path.exists(CLUST_PATH):
    with fits.open(CLUST_PATH) as hdul:
        c_d = hdul[1].data
        clust_z_map = dict(zip(c_d["TARGETID"].astype(np.int64), c_d["Z"]))

ra_l, dec_l, z_l, tid_l = [], [], [], []
for i in range(len(full_tid)):
    tid = full_tid[i]
    zv = clust_z_map.get(tid, full_data[z_col][i])
    rv, dv = full_data["RA"][i], full_data["DEC"][i]
    if np.isfinite(zv) and np.isfinite(rv) and np.isfinite(dv):
        ra_l.append(rv); dec_l.append(dv); z_l.append(zv); tid_l.append(tid)

ra, dec, z, targetid = np.array(ra_l), np.array(dec_l), np.array(z_l), np.array(tid_l)
cosmo = build_cosmology()
coords = radec_z_to_xyz(ra, dec, z, cosmo).T
tree = cKDTree(coords)
id_to_index = {int(tid): i for i, tid in enumerate(targetid)}
print(f"Success! Loaded {len(targetid)} clean objects.")

# ============================================================
# 3. SEARCH & SCAN LOGIC
# ============================================================
def build_fixed_cube(center):
    u = BASE_AXIS
    temp = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    v = np.cross(u, temp); v /= np.linalg.norm(v)
    w = np.cross(u, v)
    e1, e2, e3 = EDGE_MPC * u, EDGE_MPC * v, EDGE_MPC * w
    return np.array([center + a*e1 + b*e2 + c*e3 for a in [0,1] for b in [0,1] for c in [0,1]])

def analyze_corner(vtx):
    dist, idx = tree.query(vtx)
    label = "GALAXY" if dist < DENSITY_RADIUS_MPC else "VOID"
    return {"tid": int(targetid[idx]), "d": dist, "lab": label}

def run_search(idx):
    center = coords[idx]
    vertices = build_fixed_cube(center)
    res = [analyze_corner(v) for v in vertices]
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, f"TARGET: {targetid[idx]} | RA: {ra[idx]:.3f} DEC: {dec[idx]:.3f}\n\n")
    for i, r in enumerate(res):
        text_output.insert(tk.END, f"Corner {i}: {r['lab']} | ID: {r['tid']} | Dist: {r['d']:.2f} Mpc\n")

# ============================================================
# 4. THE PRECISION MISSILE & 3D VISUALIZER
# ============================================================
import matplotlib.pyplot as plt

last_corner_data = None
last_target_id = None

def plot_structure_3d(target_id, corner_info, coords, id_to_index, radius=150.0):
    if target_id not in id_to_index: return
    t_idx = id_to_index[target_id]
    center = coords[t_idx]
    idx_neigh = tree.query_ball_point(center, r=radius)
    neigh = coords[idx_neigh]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(neigh[:,0], neigh[:,1], neigh[:,2], s=1, alpha=0.1, color="gray")

    for c in corner_info:
        cid = c["tid"]
        ctype = c["lab"]
        pos = coords[id_to_index[cid]]
        color = "blue" if ctype == "GALAXY" else "orange"
        ax.scatter(pos[0], pos[1], pos[2], s=50, color=color)
        ax.text(pos[0], pos[1], pos[2], f"{ctype}", fontsize=8)

    ax.scatter(center[0], center[1], center[2], s=100, color="red")
    ax.set_title(f"3D Map: Target {target_id} - The Lattice View")
    plt.show()

def fire_precision_missile():
    global BASE_AXIS
    import datetime
    log_filename = "San_Salvador_Discovery.txt"
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, ">>> FLORIDA PUSH: DRILLING FOR 70 PERCENT+ RESONANCE...\n")
    root.update_idletasks()

    # Start from our current 67.0% peak
    best_rate = 67.0  
    current_best_axis = BASE_AXIS.copy()
    
    # Increase trials to 150 for the 'Florida March'
    for sweep in range(150):
        # Ultra-fine nudges (0.02) to avoid overshooting the peak
        nudge = (np.random.rand(3) - 0.5) * 0.02 
        trial_axis = current_best_axis + nudge
        trial_axis /= np.linalg.norm(trial_axis)
        
        # Test against a larger 100-galaxy sample for 'Statistically Pure' results
        sample_indices = np.random.choice(len(coords), 100, replace=False)
        total_hits = 0
        
        for idx in sample_indices:
            center = coords[idx]
            u = trial_axis
            # Calculate the cube orientation
            temp = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            v = np.cross(u, temp); v /= np.linalg.norm(v)
            w = np.cross(u, v)
            e1, e2, e3 = EDGE_MPC * u, EDGE_MPC * v, EDGE_MPC * w
            vertices = np.array([center + a*e1 + b*e2 + c*e3 for a in [0,1] for b in [0,1] for c in [0,1]])
            
            for vtx in vertices:
                dist, _ = tree.query(vtx)
                if dist < DENSITY_RADIUS_MPC:
                    total_hits += 1
        
        success_rate = (total_hits / (100 * 8)) * 100
        
        if success_rate > best_rate:
            best_rate = success_rate
            current_best_axis = trial_axis
            BASE_AXIS = current_best_axis # Lock the new steering
            text_output.insert(tk.END, f"NEW FLORIDA PEAK: {best_rate:.1f}% at Sweep {sweep}\n")
            text_output.see(tk.END)
            root.update_idletasks()
        elif sweep % 15 == 0:
            text_output.insert(tk.END, f"Marching toward Florida... Sweep {sweep}/150\n")
            root.update_idletasks()

    report = (f"\nFLORIDA MISSION REPORT\nPeak Resonance: {best_rate:.1f}%\n"
              f"Master Vector: {current_best_axis.tolist()}\n")
    
    with open(log_filename, "a", encoding="utf-8") as f:
        f.write(f"\nFLORIDA PUSH - {datetime.datetime.now()}\n{report}\n")

    text_output.insert(tk.END, report)
    messagebox.showinfo("Florida Reach", f"Final Peak: {best_rate:.1f}%")

def search_trigger():
    global last_corner_data, last_target_id
    ra_in, dec_in, z_in = entry_ra.get().strip(), entry_dec.get().strip(), entry_z.get().strip()
    try:
        target_xyz = radec_z_to_xyz(float(ra_in), float(dec_in), float(z_in), cosmo)
        _, idx = tree.query(target_xyz)
        
        center = coords[idx]
        vertices = build_fixed_cube(center)
        res = [analyze_corner(v) for v in vertices]
        last_corner_data = res
        last_target_id = int(targetid[idx])
        
        text_output.delete("1.0", tk.END)
        text_output.insert(tk.END, f"TARGET LOCKED: {last_target_id}\n\n")
        for i, r in enumerate(res):
            text_output.insert(tk.END, f"Corner {i}: {r['lab']} | Dist: {r['d']:.2f} Mpc\n")
    except: messagebox.showerror("Error", "Check RA/DEC/Z values.")

def open_3d_map():
    if last_target_id: plot_structure_3d(last_target_id, last_corner_data, coords, id_to_index)
    else: messagebox.showwarning("No Data", "Lock a target first!")

# ============================================================
# 5. THE FINAL DASHBOARD
# ============================================================
root = tk.Tk(); root.title("Cristoforo Colombo Explorer")

tk.Label(root, text="RA / DEC / Z:").grid(row=0, column=0, pady=10)
frame_gps = tk.Frame(root); frame_gps.grid(row=0, column=1)
entry_ra = tk.Entry(frame_gps, width=8); entry_ra.pack(side="left")
entry_dec = tk.Entry(frame_gps, width=8); entry_dec.pack(side="left")
entry_z = tk.Entry(frame_gps, width=8); entry_z.pack(side="left")

tk.Button(root, text="🎯 LOCK TARGET", command=search_trigger, bg="#c1f0c1", width=20).grid(row=1, column=0, columnspan=2)
text_output = tk.Text(root, width=60, height=15); text_output.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

tk.Button(root, text="🚀 FIRE PRECISION MISSILE (100 SWEEPS)", bg="red", fg="white", font=('Arial', 10, 'bold'),
          command=lambda: threading.Thread(target=fire_precision_missile, daemon=True).start()).grid(row=3, column=0, columnspan=2, pady=5)

tk.Button(root, text="🌌 SHOW 3D STRUCTURE", bg="#3498db", fg="white", font=('Arial', 10, 'bold'),
          command=open_3d_map).grid(row=4, column=0, columnspan=2, pady=5)

root.mainloop()
root = tk.Tk(); root.title("Galaxy Periodic Oscillation Explorer")

tk.Label(root, text="TARGETID:").grid(row=0, column=0, pady=5)
entry_tid = tk.Entry(root, width=30); entry_tid.grid(row=0, column=1)

tk.Label(root, text="RA / DEC / Z:").grid(row=1, column=0)
frame_gps = tk.Frame(root); frame_gps.grid(row=1, column=1)
entry_ra = tk.Entry(frame_gps, width=8); entry_ra.pack(side="left")
entry_dec = tk.Entry(frame_gps, width=8); entry_dec.pack(side="left")
entry_z = tk.Entry(frame_gps, width=8); entry_z.pack(side="left")

tk.Button(root, text="SEARCH SINGLE ID/GPS", command=search_trigger, bg="#c1f0c1").grid(row=2, column=0, columnspan=2, pady=10)

text_output = tk.Text(root, width=75, height=18); text_output.grid(row=3, column=0, columnspan=2, padx=10)

# The Precision Missile Button
btn_precision = tk.Button(root, text="🎯 FIRE PRECISION MISSILE", bg="red", fg="white", font=('Arial', 10, 'bold'),
                         command=lambda: threading.Thread(target=fire_precision_missile, daemon=True).start())
btn_precision.grid(row=5, column=0, columnspan=2, pady=5)

root.mainloop()