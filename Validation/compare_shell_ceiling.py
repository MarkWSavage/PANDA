import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# Surrounding-volume ceiling check
# =========================
# DetectorConstruction.cc auto-grows the surrounding volume as:
#   surroundingXY        = max(100 um, 1.2 x max(sensitiveXY, deadXY))
#   surroundingThickness  = max(100 um, 1.2 x (deadThickness + sensitiveThickness))
# Below ~83 um in either axis the fixed 100 um default dominates (the
# McNulty-validated regime: ~11.7 um device inside a 100 um shell, ~8.5x
# margin -- and even that generous ratio left a known ~48x under-prediction
# at high deposited energy, see Documentation/PANDA_MASTER_DESIGN Known
# Limitations). Above ~83 um the shell shrinks to a flat 1.2x-stack (10%
# per side) margin no matter how big the device gets -- the hypothesis this
# script checks is that this thinner margin under-captures nearby nuclear-
# reaction recoils, i.e. that PANDA silently gets worse, not better, above
# the ~80 um threshold.
#
# Self-consistency method (no real device exists at this size to compare
# against): at each nominal device size, PANDA was run once with the
# auto-grown shell and once with an explicit 4x-device-size reference
# shell (see Macros/run_ceiling_*.mac). If the two agree, the auto-grown
# shell is adequate at that size. If they diverge, the auto-grown shell is
# truncating real physics -- and since a 4x-wider shell is not proof of
# full convergence, only evidence of truncation if it disagrees with auto,
# a genuine gap found here would call for a convergence sweep (1.2x/2x/4x/
# 8x at one size) as a follow-up to find where the answer actually
# stabilizes.
#
# First pass (50/100/200 um, 50 MeV) found NO confirmation: all ratios
# stayed within ~1-4% of 1.0 with no growing-deficit-with-size pattern.
# This follow-up extends to 500 um-1 mm and a higher proton energy
# (200 MeV) in case the effect only appears at larger scale or where
# secondary-particle ranges are longer.
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)

qcrit_trials_fC = [50, 100, 150, 250]
nParticles = 10_000_000  # matches /run/beamOn in every run_ceiling_*.mac

# (axis_label, nominal_size_um)
SIZES = [50, 100, 200, 500, 1000]
AXES = ["XY", "T"]

# 500/1000 um points are a follow-up at a higher proton energy (200 MeV,
# matching the recoil-LET investigation's reference energy) -- the first
# pass (50/100/200 um) used 50 MeV throughout. See Macros/run_ceiling_*.mac.
ENERGY_MeV = {50: 50, 100: 50, 200: 50, 500: 200, 1000: 200}


def events_path(axis, size, condition):
    return os.path.join(
        project_root, "Results", f"Ceiling_{axis}{size}_{condition}", "events.csv"
    )


def load_probs(axis, size, condition):
    # Returns, per Qcrit threshold, both the weighted probability (used for
    # the ratio) and the raw (unweighted) event count above threshold (used
    # to estimate the ratio's statistical uncertainty via Poisson counting
    # stats -- sqrt(1/n1 + 1/n2) on the ratio of two counts). The bias
    # factor boosts real interaction rate, not the sample size, so raw
    # count -- not the weighted probability -- is what sets the noise floor.
    path = events_path(axis, size, condition)
    if not os.path.exists(path):
        print(f"[skip] {path} not found")
        return None
    df = pd.read_csv(path, usecols=["CollectedCharge_fC", "EventWeight"])
    q_fC = df["CollectedCharge_fC"].values
    w = df["EventWeight"].values
    return {
        q: {"prob": w[q_fC >= q].sum() / nParticles, "n": int((q_fC >= q).sum())}
        for q in qcrit_trials_fC
    }


# results[axis][size][condition][qcrit] = probability
results = {axis: {} for axis in AXES}
for axis in AXES:
    for size in SIZES:
        results[axis][size] = {
            "auto": load_probs(axis, size, "auto"),
            "wide": load_probs(axis, size, "wide"),
        }

# =========================
# Table: auto/wide ratio per axis x size x threshold
# =========================
print("===== AUTO-SHELL / WIDE-SHELL RATIO (1.0 = shell size doesn't matter) =====")
print("(below ~0.83x nominal-size threshold factor, both shells nearly the same size --")
print(" expect ratio near 1.0 there; departure above the ~80 um threshold confirms the ceiling)\n")

table_rows = []
for axis in AXES:
    for size in SIZES:
        auto = results[axis][size]["auto"]
        wide = results[axis][size]["wide"]
        if auto is None or wide is None:
            continue
        print(f"{axis}{size} (nominal {size} um, {ENERGY_MeV[size]} MeV):")
        for q in qcrit_trials_fC:
            pa, pw = auto[q]["prob"], wide[q]["prob"]
            na, nw = auto[q]["n"], wide[q]["n"]
            ratio = (pa / pw) if pw > 0 else float("nan")
            # Poisson counting uncertainty on the ratio of two raw counts
            # (not the bias-weighted probabilities -- the bias factor
            # boosts real interaction rate, not sample size).
            ratio_sigma = ratio * np.sqrt(1.0 / na + 1.0 / nw) if na > 0 and nw > 0 else float("nan")
            sigmas_from_1 = abs(ratio - 1.0) / ratio_sigma if ratio_sigma > 0 else float("nan")
            print(
                f"  Qcrit={q:>3} fC:  P(auto)={pa:.3e} (n={na:>7})  "
                f"P(wide)={pw:.3e} (n={nw:>7})  "
                f"ratio={ratio:.3f} +/- {ratio_sigma:.3f}  ({sigmas_from_1:.1f} sigma from 1.0)"
            )
            table_rows.append(
                {
                    "axis": axis, "size": size, "energy_MeV": ENERGY_MeV[size], "qcrit": q,
                    "P_auto": pa, "P_wide": pw, "n_auto": na, "n_wide": nw,
                    "ratio": ratio, "ratio_sigma": ratio_sigma, "sigmas_from_1": sigmas_from_1,
                }
            )
        print()

table = pd.DataFrame(table_rows)
table.to_csv(os.path.join(script_dir, "shell_ceiling_results.csv"), index=False)
print(f"Saved {os.path.join(script_dir, 'shell_ceiling_results.csv')}")

# =========================
# Verdict
# =========================
print("\n===== VERDICT =====")
print("(a departure is only meaningful if |ratio-1| clears several sigma --")
print(" small departures well within 1-2 sigma are statistical noise, not a real effect)\n")
any_significant = False
for axis in AXES:
    for size in SIZES:
        sub = table[(table["axis"] == axis) & (table["size"] == size)]
        if sub.empty:
            continue
        max_sigma_row = sub.loc[sub["sigmas_from_1"].idxmax()]
        flag = "  <-- >=3 sigma" if max_sigma_row["sigmas_from_1"] >= 3 else ""
        if flag:
            any_significant = True
        print(
            f"{axis}{size}: worst-case departure at Qcrit={int(max_sigma_row['qcrit'])} fC, "
            f"ratio={max_sigma_row['ratio']:.3f} +/- {max_sigma_row['ratio_sigma']:.3f} "
            f"({max_sigma_row['sigmas_from_1']:.1f} sigma){flag}"
        )
print(
    "\nNo threshold at any size on either axis reached 3-sigma significance."
    if not any_significant else
    "\nAt least one (axis, size, Qcrit) point reached >=3 sigma -- inspect above."
)

# =========================
# Plot: small multiples, one panel per axis x size
# =========================
fig, axes_grid = plt.subplots(2, len(SIZES), figsize=(5 * len(SIZES), 9), sharex=True, sharey=True)

for row, axis in enumerate(AXES):
    for col, size in enumerate(SIZES):
        ax = axes_grid[row][col]
        auto = results[axis][size]["auto"]
        wide = results[axis][size]["wide"]
        if auto is not None:
            ax.semilogy(qcrit_trials_fC, [auto[q]["prob"] for q in qcrit_trials_fC], 'o-', label="auto shell")
        if wide is not None:
            ax.semilogy(qcrit_trials_fC, [wide[q]["prob"] for q in qcrit_trials_fC], 's--', label="wide shell (4x)")
        ax.set_title(f"{axis} axis, nominal {size} um\n{ENERGY_MeV[size]} MeV", fontsize=10)
        ax.grid(True, which="both", alpha=0.3)

for ax in axes_grid[-1]:
    ax.set_xlabel("Qcrit (fC)")
for ax in axes_grid[:, 0]:
    ax.set_ylabel("P(collected charge >= Qcrit)")
axes_grid[0][0].legend(fontsize=8)

fig.suptitle(
    "Surrounding-volume ceiling check: auto-grown vs. 4x-wide reference shell\n"
    "convergence (overlapping curves) at 50 um expected; divergence above ~80 um "
    "would confirm the auto-grow ceiling",
    fontsize=11,
)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "PANDA_shell_ceiling.png"), dpi=300)
print(f"\nSaved {os.path.join(script_dir, 'PANDA_shell_ceiling.png')}")
