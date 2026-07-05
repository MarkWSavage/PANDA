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
E_threshold = 0.25  # MeV

# -----------------------------------
# McNulty digitized neutron data
# -----------------------------------
mcnulty_energy = np.array([
    3.2997784491440076,
    3.621027190332326,
    4.088499496475327,
    4.628197381671702,
    5.132890231621349,
    5.585740181268884,
    6.091762336354483,
    6.534420946626385,
    7.038670694864049,
    7.533172205438065,
    8.07730110775428,
    8.546988922457201,
    9.26614300100705,
    9.492124874118831,
    10.012326283987916,
    10.64773413897281
])

mcnulty_rate = np.array([
    5.825184264312032e-15,
    3.781123354537916e-15,
    1.8770107822223726e-15,
    1.218006334176351e-15,
    8.90383368650301e-16,
    6.413121772735442e-16,
    4.902219455766502e-16,
    3.860862100034026e-16,
    2.78065126745702e-16,
    2.2226603990452306e-16,
    1.673818868430848e-16,
    1.3784171397510957e-16,
    8.810335507092169e-17,
    7.366707290514941e-17,
    5.888231744032774e-17,
    4.636210454976174e-17
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
plt.title("PANDA Quenching Validation vs McNulty")

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
