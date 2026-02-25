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

def run_bulk_scan():
    import datetime
    # 1. Start the Log with the official Crew Titles
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_filename = "San_Salvador_Discovery.txt"
    
    with open(log_filename, "a") as f:
        f.write("\n" + "="*60 + "\n")
        f.write("      OFFICIAL DISCOVERY LOG: THE PERIODIC GRID\n")
        f.write(f"      Timestamp: {timestamp}\n")
        f.write("      Captain: [The Driver]\n")
        f.write("      High Navigator: Gemini\n")
        f.write("      Quartermaster: Copilot\n")
        f.write("="*60 + "\n")

    text_output.insert(tk.END, f"\n>>> NAVIGATOR: LOGGING TO {log_filename}...\n")
    root.update_idletasks() 
    
    sample_size = 100
    total_potential_hits = sample_size * 8
    total_actual_hits = 0
    
    sample_indices = np.random.choice(len(coords), sample_size, replace=False)
    
    for i, idx in enumerate(sample_indices):
        center = coords[idx]
        vertices = build_fixed_cube(center)
        cube_hits = 0
        for vtx in vertices:
            dist, _ = tree.query(vtx)
            if dist < DENSITY_RADIUS_MPC:
                total_actual_hits += 1
                cube_hits += 1
        
        # Log High-Resonance Anchors
        if cube_hits >= 5:
            with open(log_filename, "a") as f:
                f.write(f"High-Resonance Node: {targetid[idx]} | Hits: {cube_hits}/8\n")

        if i % 10 == 0:
            text_output.insert(tk.END, f"Calculating... {i}% of the haul checked\n")
            root.update_idletasks()

    avg_hits = total_actual_hits / float(sample_size)
    # The Correct Success Rate Math:
    success_rate = (total_actual_hits / float(total_potential_hits)) * 100
    
    report = (f"\n*** FINAL NAVIGATION REPORT ***\n"
              f"Average Hits: {avg_hits:.2f} per 8-corner Cube\n"
              f"Grid Success Rate: {success_rate:.1f}%\n")
    
    with open(log_filename, "a") as f:
        f.write(report + "="*60 + "\n")

    text_output.insert(tk.END, report)
    text_output.see(tk.END)
    messagebox.showinfo("Cuba Mission", "Data logged. The grid is holding steady.")
    
    text_output.insert(tk.END, report)
    text_output.see(tk.END)
    messagebox.showinfo("Cuba Mission", "Data logged. The grid is holding steady.")

# ============================================================
# 4. THE PRECISION MISSILE (HUNTING FOR 70%)
# ============================================================
def fire_precision_missile():
    import datetime
    log_filename = "San_Salvador_Discovery.txt"
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, ">>> PRECISION STRIKE: HUNTING FOR 70% RESONANCE...\n")
    root.update_idletasks()

    best_rate = 63.1  
    current_best_axis = BASE_AXIS.copy()
    
    for sweep in range(30):
        nudge = (np.random.rand(3) - 0.5) * 0.05 
        trial_axis = current_best_axis + nudge
        trial_axis /= np.linalg.norm(trial_axis)
        
        sample_indices = np.random.choice(len(coords), 50, replace=False)
        total_hits = 0
        
        for idx in sample_indices:
            center = coords[idx]
            u = trial_axis
            temp = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            v = np.cross(u, temp); v /= np.linalg.norm(v)
            w = np.cross(u, v)
            e1, e2, e3 = EDGE_MPC * u, EDGE_MPC * v, EDGE_MPC * w
            vertices = np.array([center + a*e1 + b*e2 + c*e3 for a in [0,1] for b in [0,1] for c in [0,1]])
            
            for vtx in vertices:
                dist, _ = tree.query(vtx)
                if dist < DENSITY_RADIUS_MPC:
                    total_hits += 1
        
        success_rate = (total_hits / (50 * 8)) * 100
        
        if success_rate > best_rate:
            best_rate = success_rate
            current_best_axis = trial_axis
            text_output.insert(tk.END, f"NEW PEAK FOUND: {best_rate:.1f}% Alignment!\n")
            text_output.see(tk.END)
            root.update_idletasks()
        elif sweep % 5 == 0:
            text_output.insert(tk.END, f"Sweep {sweep}/30: Adjusting orientation...\n")
            root.update_idletasks()

    report = (f"\n🏆 CUBA MISSION COMPLETE\n"
              f"Peak Grid Success: {best_rate:.1f}%\n"
              f"Final Master Axis: {current_best_axis.tolist()}\n")
    
    with open(log_filename, "a") as f:
        f.write(f"\nPRECISION SWEEP - {datetime.datetime.now()}\n")
        f.write(report + "="*60 + "\n")

    text_output.insert(tk.END, report)
    messagebox.showinfo("Mission Success", f"New Peak: {best_rate:.1f}%")
    
def search_trigger():
    tid_in = entry_tid.get().strip()
    ra_in, dec_in, z_in = entry_ra.get().strip(), entry_dec.get().strip(), entry_z.get().strip()
    
    if tid_in.isdigit() and int(tid_in) in id_to_index:
        run_search(id_to_index[int(tid_in)])
    elif ra_in and dec_in and z_in:
        try:
            target_xyz = radec_z_to_xyz(float(ra_in), float(dec_in), float(z_in), cosmo)
            _, idx = tree.query(target_xyz)
            run_search(idx)
        except: messagebox.showerror("Error", "Check RA/DEC/Z values.")
    else: messagebox.showerror("Error", "Enter ID or RA/DEC/Z")

# ============================================================
# 4. THE GUI (The Dashboard)
# ============================================================
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