import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


RESULTS_DIR = os.path.join("Results", "Current")
os.makedirs(RESULTS_DIR, exist_ok=True)


# ----------------------------
# Load PANDA output
# ----------------------------
csv_file = os.path.join(RESULTS_DIR, "events.csv")

print("Loading:", csv_file)

data = pd.read_csv(csv_file)

total_keV = data["Total_keV"].values
nonzero_keV = total_keV[total_keV > 0]

proton_keV   = data["Proton_keV"].values
electron_keV = data["Electron_keV"].values
recoil_keV   = data["Recoil_keV"].values

# ----------------------------
# Filter Data
# ----------------------------
proton_nonzero   = proton_keV[proton_keV > 0]
electron_nonzero = electron_keV[electron_keV > 0]
recoil_nonzero   = recoil_keV[recoil_keV > 0]

# ----------------------------
# Basic statistics
# ----------------------------
print("\nEnergy statistics:")
print(f"Min:  {np.min(total_keV):.3f} keV")
print(f"Mean: {np.mean(total_keV):.3f} keV")
print(f"Max:  {np.max(total_keV):.3f} keV")
print(f"Total events: {len(total_keV)}")
print(f"Triggered events: {len(nonzero_keV)}")
print(f"Efficiency: {len(nonzero_keV)/len(total_keV):.8f}")

# ----------------------------
# Histogram
# ----------------------------
bins = np.linspace(
    0,
    np.max(nonzero_keV),
    200
)

hist_total, edges = np.histogram(
    nonzero_keV,
    bins=bins
)

hist_proton, _ = np.histogram(
    proton_nonzero,
    bins=bins
)

hist_electron, _ = np.histogram(
    electron_nonzero,
    bins=bins
)

hist_recoil, _ = np.histogram(
    recoil_nonzero,
    bins=bins
)

centers = 0.5 * (edges[:-1] + edges[1:])


# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(10, 6))

plt.step(
    centers,
    hist_total,
    where="mid",
    label="Total"
)

plt.step(
    centers,
    hist_proton,
    where="mid",
    label="Proton"
)

plt.step(
    centers,
    hist_electron,
    where="mid",
    label="Electron"
)

plt.step(
    centers,
    hist_recoil,
    where="mid",
    label="Recoil"
)

plt.legend()
plt.xlabel("Deposited Energy (keV)")
plt.ylabel("Counts")
plt.title("PANDAEX Raw Deposited Energy Spectrum")
plt.yscale("log")

outfile = os.path.join(
    RESULTS_DIR,
    "PANDAEX_component_spectrum.png"
)

plt.savefig(outfile, dpi=300)

print("\nSaved:")
print(outfile)

plt.show()
