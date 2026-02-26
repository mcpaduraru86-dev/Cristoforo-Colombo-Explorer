import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.fftpack import fft, fftfreq

# --- 1. THE LOADING DOCK ---
# Make sure "full.dat.fits" is in the same folder as this script!
file_name = "full.dat.fits" 

print(f"🚛 ATTEMPTING TO LOAD: {file_name}")

try:
    # This creates the 'df' that the error was complaining about
    dat = Table.read(file_name, format='fits')
    df = dat.to_pandas()
    print(f"✅ CARGO LOADED: {len(df)} galaxies found.")
except Exception as e:
    print(f"❌ BREAKDOWN: Could not find the file. Error: {e}")
    exit()

# --- 2. THE BORE-SIGHT (V23.0) ---
print(">>> V23.0: BORING ALONG THE COSMIC AXIS (RA 168 <-> 348)...")

# Now 'df' exists, so the scan can proceed
shaft_mask = ((df['RA'] > 158) & (df['RA'] < 178)) | ((df['RA'] > 338) & (df['RA'] < 358))
shaft_z = df[shaft_mask]['Z']

# --- 3. THE MATH (Volvo Momentum) ---
# Simple distance calc for the scan
H0 = 67.4
c = 299792.458
dist = (c * shaft_z) / H0

counts, bins = np.histogram(dist, bins=500, range=(0, 600))
yf = np.abs(fft(counts - np.mean(counts)))
xf = 1/fftfreq(len(yf), (bins[1] - bins[0]))

# --- 4. THE DISPLAY ---
plt.figure(figsize=(10, 6))
plt.plot(xf[(xf > 40) & (xf < 500)], yf[(xf > 40) & (xf < 500)], color='blue')
plt.axvline(137.5, color='red', linestyle='--', label='137.5 Mpc Tesseract Beat')
plt.title("DRIVE-SHAFT SCAN: PERIODIC OSCILLATIONS")
plt.xlabel("Distance Scale (Mpc)")
plt.ylabel("Signal Strength")
plt.legend()

print("🏁 SCAN COMPLETE. Opening the display window...")
plt.show()