import os
import threading
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
from astropy.io import fits
from scipy.spatial import cKDTree

from utils import build_cosmology, radec_z_to_xyz
from plotting import plot_real_and_hypothetical

# ============================================================
# CONFIG
# ============================================================
FULL_PATH = r"C:\Users\Marian Catalin\Documents\Ianus Paduraru\BGS_BRIGHT_full.dat.fits"
CLUST_PATH = r"C:\Users\Marian Catalin\Documents\Ianus Paduraru\BGS_BRIGHT_NGC_clustering.dat.fits"

DENSITY_RADIUS_MPC = 10.0
EDGE_MPC = 137.5
BASE_AXIS = np.array([0.303, 0.808, 0.505], dtype=float)
BASE_AXIS /= np.linalg.norm(BASE_AXIS)

# ============================================================
# DATA LOADING
# ============================================================
print("Loading 5GB Catalog...")
with fits.open(FULL_PATH) as hdul:
    full_data = hdul[1].data
    full_tid = full_data["TARGETID"].astype(np.int64)
    possible_z = ['Z', 'Z_not4clus', 'Z_HP', 'Z_PHOT']
    z_col = next((n for n in possible_z if n in full_data.names), 'Z')

clust_z_map = {}
if os.path.exists(CLUST_PATH):
    with fits.open(CLUST_PATH) as hdul:
        c_d = hdul[1].data
        clust_z_map = dict(zip(c_d["TARGETID"].astype(np.int64), c_d["Z"]))

ra_l, dec_l, z_l, tid_l = [], [], [], []
for i in range(len(full_tid)):
    tid = full_tid[i]
    zv, rv, dv = clust_z_map.get(tid, full_data[z_col][i]), full_data["RA"][i], full_data["DEC"][i]
    if np.isfinite(zv) and np.isfinite(rv) and np.isfinite(dv):
        ra_l.append(rv); dec_l.append(dv); z_l.append(zv); tid_l.append(tid)

ra, dec, z, targetid = np.array(ra_l), np.array(dec_l), np.array(z_l), np.array(tid_l)
cosmo = build_cosmology()
Dc = cosmo.comoving_distance(z).value
ra_rad, dec_rad = np.deg2rad(ra), np.deg2rad(dec)
x = Dc * np.cos(dec_rad) * np.cos(ra_rad)
y = Dc * np.cos(dec_rad) * np.sin(ra_rad)
z3d = Dc * np.sin(dec_rad)
coords = np.vstack([x, y, z3d]).T
tree = cKDTree(coords)
id_to_index = {int(tid): i for i, tid in enumerate(targetid)}
print(f"Loaded {len(targetid)} galaxies.")

# ============================================================
# SEARCH ENGINE
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
    neigh = tree.query_ball_point(vtx, DENSITY_RADIUS_MPC)
    label = "GALAXY" if dist < DENSITY_RADIUS_MPC else "VOID"
    return {"v": vtx, "idx": idx, "tid": int(targetid[idx]), "d": dist, "lab": label, "c": len(neigh)}

def run_search(idx):
    center = coords[idx]
    vertices = build_fixed_cube(center)
    res = [analyze_corner(v) for v in vertices]
    text_output.delete("1.0", tk.END)
    text_output.insert(tk.END, f"TARGET: {targetid[idx]} | RA: {ra[idx]:.3f} DEC: {dec[idx]:.3f}\n\n")
    for i, r in enumerate(res):
        text_output.insert(tk.END, f"Corner {i}: {r['lab']} | ID: {r['tid']} | Dist: {r['d']:.2f} Mpc\n")
    btn_plot.config(state="normal", command=lambda: plot_real_and_hypothetical(coords, [idx], vertices, None, [r["idx"] for r in res], None, targetid))

def search_trigger():
    tid_in = entry_tid.get().strip()
    ra_in, dec_in, z_in = entry_ra.get().strip(), entry_dec.get().strip(), entry_z.get().strip()
    
    if tid_in.isdigit() and int(tid_in) in id_to_index:
        run_search(id_to_index[int(tid_in)])
    elif ra_in and dec_in and z_in:
        try:
            # RA/DEC to XYZ Search
            tx, ty, tz = radec_z_to_xyz(float(ra_in), float(dec_in), float(z_in), cosmo)
            _, idx = tree.query([tx, ty, tz])
            run_search(idx)
        except: messagebox.showerror("Error", "Check RA/DEC/Z values.")
    else: messagebox.showerror("Error", "Enter ID or RA/DEC/Z")

# ============================================================
# THE GUI (With all boxes restored)
# ============================================================
root = tk.Tk(); root.title("Galaxy Cube Explorer")

# Layout
tk.Label(root, text="TARGETID:").grid(row=0, column=0)
entry_tid = tk.Entry(root); entry_tid.grid(row=0, column=1)

tk.Label(root, text="RA / DEC / Z:").grid(row=1, column=0)
frame_gps = tk.Frame(root); frame_gps.grid(row=1, column=1)
entry_ra = tk.Entry(frame_gps, width=8); entry_ra.pack(side="left")
entry_dec = tk.Entry(frame_gps, width=8); entry_dec.pack(side="left")
entry_z = tk.Entry(frame_gps, width=8); entry_z.pack(side="left")

tk.Button(root, text="Search", command=search_trigger).grid(row=2, column=0, columnspan=2)
btn_plot = tk.Button(root, text="Plot Cube", state="disabled"); btn_plot.grid(row=3, column=0, columnspan=2)
text_output = tk.Text(root, width=70, height=15); text_output.grid(row=4, column=0, columnspan=2)

root.mainloop()