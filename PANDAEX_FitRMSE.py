import os
import pandas as pd
import numpy as np


RESULTS_DIR = os.path.join("Results", "Current")

csv_file = os.path.join(RESULTS_DIR, "events.csv")

print("Loading:", csv_file)

data = pd.read_csv(csv_file)

N_total = len(data)

proton_keV = data["Proton_keV"].values
electron_keV = data["Electron_keV"].values
recoil_keV = data["Recoil_keV"].values

# -----------------------------
# McNulty digitized data
# -----------------------------
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

threshold = 0.25

quench_scan = [0.55, 0.60, 0.65, 0.70]

print("\nRMSE(log10) results:")
print("---------------------")

for k in quench_scan:

    Edep_keV = (
        proton_keV +
        electron_keV +
        k * recoil_keV
    )

    Edep_MeV = Edep_keV / 1000.0

    Edep_MeV = Edep_MeV[Edep_MeV >= threshold]

    cumulative_prob = []

    for Eth in mcnulty_energy:
        count = np.sum(Edep_MeV >= Eth)
        cumulative_prob.append(count / N_total)

    cumulative_prob = np.array(cumulative_prob)

    # normalize at 5 MeV
    ref_E = 5.0

    panda_ref = np.interp(
        ref_E,
        mcnulty_energy,
        cumulative_prob
    )

    mcnulty_ref = np.interp(
        ref_E,
        mcnulty_energy,
        mcnulty_rate
    )

    norm_factor = mcnulty_ref / panda_ref

    panda_norm = cumulative_prob * norm_factor

    rmse = np.sqrt(
        np.mean(
            (
                np.log10(panda_norm) -
                np.log10(mcnulty_rate)
            ) ** 2
        )
    )

    print(f"Quench {k:.2f} -> RMSE = {rmse:.5f}")

