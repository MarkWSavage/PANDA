import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


RESULTS_DIR = os.path.join("Results", "Current")
csv_file = os.path.join(RESULTS_DIR, "events.csv")

data = pd.read_csv(csv_file)

# -----------------------------
# Raw charge from PANDA
# -----------------------------
Q_raw = data["Charge_fC"].values

N = len(Q_raw)

# -----------------------------
# Suppression parameter
# -----------------------------
Q0_scan = [5, 10, 25, 50, 100]

plt.figure(figsize=(10, 7))

# -----------------------------
# Raw curve
# -----------------------------
Q_thresholds = np.logspace(
    np.log10(max(1e-3, np.min(Q_raw[Q_raw > 0]))),
    np.log10(np.max(Q_raw)),
    300
)

sigma_raw = []

for Qth in Q_thresholds:
    sigma_raw.append(np.sum(Q_raw >= Qth) / N)

sigma_raw = np.array(sigma_raw)

plt.loglog(
    Q_thresholds,
    sigma_raw,
    linewidth=3,
    label="Raw PANDA"
)

# -----------------------------
# Suppressed curves
# -----------------------------
for Q0 in Q0_scan:

    Q_eff = Q_raw * (1.0 - np.exp(-Q_raw / Q0))

    sigma_eff = []

    for Qth in Q_thresholds:
        sigma_eff.append(np.sum(Q_eff >= Qth) / N)

    sigma_eff = np.array(sigma_eff)

    plt.loglog(
        Q_thresholds,
        sigma_eff,
        linewidth=2,
        label=f"Suppressed Q0={Q0} fC"
    )

plt.xlabel("Collected Charge (fC)")
plt.ylabel("Upset Probability")
plt.title("PANDA Recombination Suppression Model")
plt.grid(True, which="both")
plt.legend()

outfile = os.path.join(
    RESULTS_DIR,
    "PANDAEX_recombination_compare.png"
)

plt.savefig(outfile, dpi=300)
plt.close()

print("Saved:")
print(outfile)
