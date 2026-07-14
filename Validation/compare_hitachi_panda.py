import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# Hitachi HM68512 measured data
# =========================
# LaBel, Moran, Hawkins, Cooley, Seidleck, Gates, Smith,
# Stassinopoulos, Marshall, Dale (NASA GSFC), "Single Event Effect
# Proton and Heavy Ion Test Results for Candidate Spacecraft
# Electronics," NSREC-era report (~1994), openly available at
# https://nepp.nasa.gov/docuploads/8742E655-056B-4568-9F06CCA49C54B502/nsrec94a.pdf
#
# Hitachi HM68512, 4 Mbit SRAM, CMOS-on-EPI process. Proton-tested at
# UC Davis, 22-63 MeV, dynamic and static modes. The report states
# "little energy dependence" over that range and gives two device-
# level saturation cross sections (confirmed by direct text
# extraction, not a digitized plot):
static_cross_section_cm2 = 1e-6      # cm^2/device
dynamic_cross_section_cm2 = 2.2e-5   # cm^2/device
tested_energy_range_MeV = (22, 63)
N_BITS = 4 * 1024 * 1024  # 4 Mbit

# =========================
# PANDA data
# =========================
# Unlike the McNulty/CUPID comparison (compare_mcnulty_panda.py),
# where the sensitive-volume geometry is stated exactly in the source
# report, the HM68512's geometry and critical charge are NOT published
# anywhere found. Neither is scanned freely: geometry (sensitiveXY,
# sensitiveThickness) requires re-simulation (see Macros/run_hitachi_*
# .mac and run_hitachi_XY*_*MeV.mac), so this covers a coarse 2x2
# corner scan (small/large footprint x thin/thick collection depth)
# around the original 3x3x3 um baseline guess, at 4 energies each
# (15/22/50/80 MeV -- below, bottom-of, mid, and above the tested
# 22-63 MeV range). Critical charge (Qcrit) is swept for free from the
# recorded CollectedCharge_fC distribution, no re-simulation needed.
# deadThickness is fixed at 5 um SiO2 throughout -- not itself
# scanned, a remaining open assumption.
#
# This is NOT a clean pass/fail check -- it's a plausibility/shape
# scan across the geometry+Qcrit space we don't have measured values
# for. The primary signal is which (geometry, Qcrit) combinations
# reproduce both the reported magnitude (~1e-6 cm^2/device static) AND
# the reported flatness across 22-63 MeV, and whether that's mostly
# controlled by Qcrit (as the single-geometry run suggested) or
# whether geometry matters too.
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
nParticles = 10_000_000  # matches /run/beamOn in every run_hitachi_*.mac

qcrit_trials_fC = [50, 100, 150, 250]

# (label, sensitiveXY_um, sensitiveThickness_um, results_tag, energies)
# results_tag=None -> baseline dir naming Results/Hitachi_<E>MeV
# results_tag=str  -> Results/Hitachi_<tag>_<E>MeV
GEOMETRIES = [
    ("XY3T3 (baseline)", 3, 3, None,     [15, 22, 30, 40, 50, 63, 80]),
    ("XY2T1 (small/thin)", 2, 1, "XY2T1", [15, 22, 50, 80]),
    ("XY2T5 (small/thick)", 2, 5, "XY2T5", [15, 22, 50, 80]),
    ("XY5T1 (large/thin)", 5, 1, "XY5T1", [15, 22, 50, 80]),
    ("XY5T5 (large/thick)", 5, 5, "XY5T5", [15, 22, 50, 80]),
]


def events_path(tag, E):
    dirname = f"Hitachi_{E}MeV" if tag is None else f"Hitachi_{tag}_{E}MeV"
    return os.path.join(project_root, "Results", dirname, "events.csv")


# geometry_results[label] = {qcrit: {E: sigma_device}}
geometry_results = {}

for label, xy, t, tag, energies in GEOMETRIES:
    beam_area_cm2 = (xy * 1e-4) ** 2
    per_qcrit = {q: {} for q in qcrit_trials_fC}

    for E in energies:
        path = events_path(tag, E)
        if not os.path.exists(path):
            print(f"[skip] {path} not found")
            continue
        df = pd.read_csv(path, usecols=["CollectedCharge_fC", "EventWeight"])
        q_fC = df["CollectedCharge_fC"].values
        w = df["EventWeight"].values
        for q in qcrit_trials_fC:
            prob = w[q_fC >= q].sum() / nParticles
            sigma_bit = beam_area_cm2 * prob
            per_qcrit[q][E] = N_BITS * sigma_bit

    geometry_results[label] = per_qcrit

# =========================
# Best-fit Qcrit per geometry (min log-magnitude error vs static)
# =========================
print("===== BEST-FIT Qcrit PER GEOMETRY (vs measured static, 1e-6 cm^2/device) =====")
best_qcrit = {}
for label, xy, t, tag, energies in GEOMETRIES:
    per_qcrit = geometry_results[label]
    best_q, best_err = None, np.inf
    for q in qcrit_trials_fC:
        vals = np.array([per_qcrit[q][E] for E in energies if E in per_qcrit[q]])
        if vals.size == 0 or np.any(vals <= 0):
            continue
        log_err = np.sqrt(np.mean((np.log10(vals) - np.log10(static_cross_section_cm2)) ** 2))
        if log_err < best_err:
            best_err, best_q = log_err, q
    best_qcrit[label] = best_q

    print(f"\n{label} (sensitiveXY={xy} um, sensitiveThickness={t} um):")
    if best_q is None:
        print("  no valid data")
        continue
    vals = np.array([per_qcrit[best_q][E] for E in energies if E in per_qcrit[best_q]])
    ratio = vals / static_cross_section_cm2
    flatness = vals.max() / vals.min() if vals.size > 1 else float("nan")
    print(f"  best Qcrit = {best_q} fC (log RMSE = {best_err:.2f} decades)")
    print(f"  PANDA/measured-static ratio: {ratio.min():.2f}x to {ratio.max():.2f}x")
    print(f"  PANDA max/min across available energies (flatness, "
          f"measured is ~1x over 22-63 MeV): {flatness:.2f}x")

# =========================
# Plot: small multiples, one panel per geometry
# =========================
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
axes = axes.flatten()

for ax, (label, xy, t, tag, energies) in zip(axes, GEOMETRIES):
    per_qcrit = geometry_results[label]
    for q in qcrit_trials_fC:
        es = [E for E in energies if E in per_qcrit[q]]
        vals = [per_qcrit[q][E] for E in es]
        if es:
            style = '-' if q != best_qcrit[label] else '-'
            lw = 3 if q == best_qcrit[label] else 1.5
            ax.semilogy(es, vals, 'o-', label=f"Qcrit={q} fC", linewidth=lw)

    ax.axhline(static_cross_section_cm2, color="black", linestyle="-", linewidth=1)
    ax.axhline(dynamic_cross_section_cm2, color="black", linestyle="--", linewidth=1)
    ax.axvspan(*tested_energy_range_MeV, color="gray", alpha=0.15)
    ax.set_title(f"{label}\nXY={xy} um, T={t} um", fontsize=10)
    ax.grid(True, which="both", alpha=0.3)

axes[0].legend(fontsize=7, loc="lower right")
for ax in axes[3:]:
    ax.set_xlabel("Proton Energy (MeV)")
for ax in [axes[0], axes[3]]:
    ax.set_ylabel("Cross Section (cm$^2$/device)")
axes[-1].axis("off")

fig.suptitle(
    "PANDA vs. Hitachi HM68512 (LaBel et al. 1994) -- geometry x Qcrit scan\n"
    "solid black = measured static (1e-6), dashed black = measured dynamic (2.2e-5), "
    "shaded = tested 22-63 MeV",
    fontsize=11,
)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "PANDA_vs_Hitachi.png"), dpi=300)
print(f"\nSaved {os.path.join(script_dir, 'PANDA_vs_Hitachi.png')}")
