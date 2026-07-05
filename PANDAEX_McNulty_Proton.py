import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


RESULTS_DIR = os.path.join("Results", "Current")
os.makedirs(RESULTS_DIR, exist_ok=True)

# -----------------------------------
# Load PANDA data
# -----------------------------------
csv_file = os.path.join(RESULTS_DIR, "events.csv")

print("Loading:", csv_file)

data = pd.read_csv(csv_file)

N_total = len(data)

print(f"\nTotal primaries: {N_total}")

# -----------------------------------
# Load PANDA components
# -----------------------------------
proton_keV = data["Proton_keV"].values
electron_keV = data["Electron_keV"].values
recoil_keV = data["Recoil_keV"].values

# -----------------------------------
# Quenching scan
# -----------------------------------
quench_scan = [1.0, 0.8, 0.6, 0.4, 0.2]

# Fixed detector threshold
E_threshold = 0.40  # MeV

# -----------------------------------
# McNulty digitized neutron data
# -----------------------------------
mcnulty_energy = np.array([
    0.05007049345417869,
    2.0524471299093654,
    3.116334340382678,
    5.137764350453173,
    7.049748237663646,
    8.936032225579053
])

mcnulty_rate = np.array([
    9.421418903256075e-15,
    1.55099994522946e-15,
    5.9730115353128e-16,
    2.8720689646078423e-16,
    7.174273775007409e-17,
    1.1637940511694179e-17
])

# -----------------------------------
# Energy thresholds for cumulative plot
# -----------------------------------
thresholds = np.logspace(
    np.log10(E_threshold),
    np.log10(15.0),
    300
)

# -----------------------------------
# Reference normalization point
# -----------------------------------
ref_E = 5.0

mcnulty_ref = np.interp(
    ref_E,
    mcnulty_energy,
    mcnulty_rate
)

# -----------------------------------
# Plot
# -----------------------------------
plt.figure(figsize=(8, 6))

for k in quench_scan:

    # Effective ionizing energy
    Edep_keV = (
        proton_keV +
        electron_keV +
        k * recoil_keV
    )

    Edep_MeV = Edep_keV / 1000.0

    # Remove zeroes
    Edep_MeV = Edep_MeV[Edep_MeV > 0]

    # Apply threshold
    Edep_MeV = Edep_MeV[Edep_MeV >= E_threshold]

    cumulative_prob = []

    for Eth in thresholds:
        count = np.sum(Edep_MeV >= Eth)
        cumulative_prob.append(count / N_total)

    cumulative_prob = np.array(cumulative_prob)

    # Normalize to McNulty at 5 MeV
    panda_ref = np.interp(
        ref_E,
        thresholds,
        cumulative_prob
    )

    norm_factor = mcnulty_ref / panda_ref

    cumulative_prob_norm = cumulative_prob * norm_factor

    print(f"Quench = {k:.2f}, norm factor = {norm_factor:.3e}")

    plt.loglog(
        thresholds,
        cumulative_prob_norm,
        linewidth=2,
        label=f"PANDA (quench={k})"
    )

# Plot McNulty
plt.loglog(
    mcnulty_energy,
    mcnulty_rate,
    'x-',
    linewidth=2,
    label="McNulty (1986)"
)

plt.xlabel("Energy Deposited (MeV)")
plt.ylabel("Burst Generation Rate")
plt.title("PANDA Proton Validation vs McNulty (37 MeV)")

plt.grid(True, which="both")
plt.legend()

outfile = os.path.join(
    RESULTS_DIR,
    "PANDAEX_McNulty_quench.png"
)

plt.savefig(outfile, dpi=300)

print("\nSaved:")
print(outfile)

plt.show()
