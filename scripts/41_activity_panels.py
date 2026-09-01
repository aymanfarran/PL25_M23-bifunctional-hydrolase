#!/usr/bin/env python3
"""Figures 6 and 7 - wet-lab validation of PL25_M23 lytic activity.

Figure 6  killing activity
    a  primary antibacterial assay   (100 ug/mL, 60 min)
    c  time-kill kinetics            (100 ug/mL, 0-180 min)
    d  dose-response                 (0-300 ug/mL, 60 min)
    Panel b of the published figure is a photograph of the drop plates and is
    composited in by hand; it is not produced here.

Figure 7  biochemical properties
    a  metal dependence   (1 mM EDTA +/- 100 uM ZnCl2, 100 ug/mL, 60 min)
    b  salt tolerance     (0 / 0.5 / 1.0 M NaCl, 100 ug/mL, 60 min)

Monochrome throughout: strains are distinguished by marker shape in the line
plots and by fill/hatch in the grouped bars, so both figures survive greyscale
reproduction. No figure-level title or footnote is drawn - assay conditions,
replicate structure and significance thresholds belong to the manuscript
figure captions.

Statistics: Welch's t-test on log10-transformed CFU/mL (untransformed counts
are strongly right-skewed), n = 3 technical replicates, single biological
experiment.

Inputs:
    results/20_activity/turbidity_raw.tsv
    results/20_activity/dose_response_raw.tsv
    results/20_activity/time_kill_raw.tsv
    results/20_activity/edta_raw.tsv
    results/20_activity/salt_raw.tsv

Outputs:
    figures/06_activity_killing.{pdf,svg,png}
    figures/07_activity_biochemistry.{pdf,svg,png}
    results/20_activity/figure_stats_summary.tsv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "20_activity"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

INK = "black"
GREY_UNTREATED = "#BFBFBF"
GREY_GRID = "#E8E8E8"

STRAINS = ["Bacillus subtilis", "Staphylococcus warneri",
           "MRSA", "Virgibacillus salarius PL25"]
MARKER = {"Bacillus subtilis": "o", "Staphylococcus warneri": "s",
          "MRSA": "D", "Virgibacillus salarius PL25": "^"}
HATCH = {"Bacillus subtilis": "", "Staphylococcus warneri": "///",
         "MRSA": "xxx", "Virgibacillus salarius PL25": "..."}
FILL = {"Bacillus subtilis": "white", "Staphylococcus warneri": "#F0F0F0",
        "MRSA": "#B8B8B8", "Virgibacillus salarius PL25": "#606060"}
LABEL = {"Bacillus subtilis": r"$\it{B.\ subtilis}$",
         "Staphylococcus warneri": r"$\it{S.\ warneri}$",
         "MRSA": "MRSA",
         "Virgibacillus salarius PL25": r"$\it{V.\ salarius}$ PL25"}

REPS = ["rep1_CFU_per_mL", "rep2_CFU_per_mL", "rep3_CFU_per_mL"]

stat_rows = []   # accumulates every test for the audit table


def log_cfu(df_row):
    """log10 of the three replicate CFU/mL values in one row."""
    return np.log10(df_row[REPS].astype(float).values)


def welch(a, b):
    """Welch's t-test on two log10-CFU vectors; returns the p-value."""
    try:
        return float(stats.ttest_ind(a, b, equal_var=False).pvalue)
    except Exception:
        return np.nan


def stars(p):
    """Significance symbol for a p-value, following the usual convention."""
    if p is None or np.isnan(p):
        return ""
    for cutoff, symbol in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")):
        if p <= cutoff:
            return symbol
    return "ns"


def record(figure, panel, strain, comparison, reduction, p):
    stat_rows.append({"figure": figure, "panel": panel, "strain": strain,
                      "comparison": comparison, "log10_reduction": reduction,
                      "p_welch": p, "symbol": stars(p), "n": 3})


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", ls="-", lw=0.5, color=GREY_GRID, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=9, colors=INK)


def panel_letter(ax, letter):
    ax.text(-0.11, 1.06, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", color=INK, ha="left", va="bottom")


def threshold_line(ax, x_text, label="3-log"):
    """Dashed marker at the 3-log10 (99.9 %) reduction level."""
    ax.axhline(3, color="#888", ls=(0, (3, 3)), lw=0.7, zorder=1)
    ax.text(x_text, 3.05, label, fontsize=8, color=INK, va="bottom")


plt.rcParams.update({
    "font.family": "DejaVu Sans", "text.color": INK,
    "axes.edgecolor": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
    "axes.linewidth": 0.9, "hatch.linewidth": 0.6,
})

# ── Load ─────────────────────────────────────────────────────────────
turbidity = pd.read_csv(DATA / "turbidity_raw.tsv", sep="\t")
dose = pd.read_csv(DATA / "dose_response_raw.tsv", sep="\t")
time_kill = pd.read_csv(DATA / "time_kill_raw.tsv", sep="\t")
edta = pd.read_csv(DATA / "edta_raw.tsv", sep="\t")
salt = pd.read_csv(DATA / "salt_raw.tsv", sep="\t")

DOSES = [0, 50, 100, 300]
TIMEPOINTS = [0, 30, 60, 120, 180]
NACL = [0.0, 0.5, 1.0]


# ═════════════════════════════════════════════════════════════════════
# FIGURE 6 - killing activity
# ═════════════════════════════════════════════════════════════════════
fig6 = plt.figure(figsize=(13.5, 8.5))
gs6 = fig6.add_gridspec(2, 2, hspace=0.55, wspace=0.30)

# ── a  primary antibacterial assay ───────────────────────────────────
ax = fig6.add_subplot(gs6[0, 0])
width = 0.36
for i, strain in enumerate(STRAINS):
    sub = turbidity[turbidity.strain == strain]
    control = log_cfu(sub[sub.treatment == "Untreated"].iloc[0])
    treated = log_cfu(sub[sub.treatment == "PL25_M23_100ugmL"].iloc[0])
    p = welch(control, treated)
    reduction = control.mean() - treated.mean()
    record("6", "a", strain, "treated vs untreated", reduction, p)

    ax.bar(i - width / 2, control.mean(), width, yerr=control.std(ddof=1),
           capsize=3, color=GREY_UNTREATED, edgecolor=INK, lw=0.6, ecolor=INK,
           zorder=3, label="Untreated" if i == 0 else None)
    ax.bar(i + width / 2, treated.mean(), width, yerr=treated.std(ddof=1),
           capsize=3, color=INK, edgecolor=INK, lw=0.6, ecolor=INK,
           zorder=3, label="+ PL25_M23" if i == 0 else None)

    # Significance bracket spanning the pair, with the reduction above it.
    top = max(control.mean() + control.std(ddof=1),
              treated.mean() + treated.std(ddof=1))
    bracket = top + 0.35
    ax.plot([i - width / 2, i + width / 2], [bracket, bracket], color=INK, lw=0.8)
    symbol = stars(p)
    ax.text(i, bracket + 0.05, symbol, ha="center", va="bottom",
            fontsize=11 if symbol != "ns" else 8.5,
            fontweight="bold" if symbol != "ns" else "normal", color=INK)
    ax.text(i, bracket + 1.05, f"Δ {reduction:.2f}",
            ha="center", fontsize=8.5, color=INK)

ax.set_xticks(range(len(STRAINS)))
ax.set_xticklabels([LABEL[s] for s in STRAINS], fontsize=9)
ax.set_ylabel("log$_{10}$ CFU/mL", fontsize=10)
ax.set_ylim(0, 11.5)
ax.legend(loc="upper right", fontsize=8.5, frameon=False)
style(ax)
panel_letter(ax, "a")

# Panel b of the published figure is the drop-plate photograph, composited
# separately; the slot is left empty here.
fig6.add_subplot(gs6[0, 1]).axis("off")

# ── c  time-kill kinetics ────────────────────────────────────────────
ax = fig6.add_subplot(gs6[1, 0])
for strain in STRAINS:
    sub = time_kill[time_kill.strain == strain]
    means, sds = [], []
    for t in TIMEPOINTS:
        control = log_cfu(sub[(sub.treatment == "Untreated") & (sub.time_min == t)].iloc[0])
        treated = log_cfu(sub[(sub.treatment == "PL25_M23") & (sub.time_min == t)].iloc[0])
        reductions = control - treated
        p = welch(control, treated)
        record("6", "c", strain, f"{t} min vs untreated", reductions.mean(), p)
        means.append(reductions.mean())
        sds.append(reductions.std(ddof=1))
        symbol = stars(p)
        if symbol not in ("", "ns") and reductions.mean() > 0.5:
            ax.text(t, reductions.mean() + reductions.std(ddof=1) + 0.1, symbol,
                    ha="center", va="bottom", fontsize=9, fontweight="bold", color=INK)
    ax.errorbar(TIMEPOINTS, means, yerr=sds, marker=MARKER[strain], ms=9,
                capsize=3.5, color=INK, mfc=INK, mec=INK, mew=0.9, lw=1.6,
                ecolor=INK, label=LABEL[strain], zorder=4)

threshold_line(ax, 182)
ax.set_xlabel("Time (min)", fontsize=10)
ax.set_ylabel("log$_{10}$ CFU reduction", fontsize=10)
ax.set_xticks(TIMEPOINTS)
ax.set_ylim(-0.3, 5.5)
ax.legend(loc="lower right", fontsize=8.5, frameon=False, ncol=2,
          handletextpad=0.4, columnspacing=0.8, labelcolor=INK)
style(ax)
panel_letter(ax, "c")

# ── d  dose-response ─────────────────────────────────────────────────
ax = fig6.add_subplot(gs6[1, 1])
for strain in STRAINS:
    sub = dose[dose.strain == strain]
    means, sds = [], []
    for d in DOSES:
        control = log_cfu(sub[(sub.treatment == "Untreated") & (sub.dose_ug_mL == d)].iloc[0])
        treated = log_cfu(sub[(sub.treatment == "PL25_M23") & (sub.dose_ug_mL == d)].iloc[0])
        reductions = control - treated
        p = welch(control, treated)
        record("6", "d", strain, f"{d} ug/mL vs untreated", reductions.mean(), p)
        means.append(reductions.mean())
        sds.append(reductions.std(ddof=1))
        symbol = stars(p)
        if symbol not in ("", "ns") and reductions.mean() > 0.5:
            ax.text(d, reductions.mean() + reductions.std(ddof=1) + 0.1, symbol,
                    ha="center", va="bottom", fontsize=9, fontweight="bold", color=INK)
    ax.errorbar(DOSES, means, yerr=sds, marker=MARKER[strain], ms=9,
                capsize=3.5, color=INK, mfc=INK, mec=INK, mew=0.9, lw=1.6,
                ecolor=INK, label=LABEL[strain], zorder=4)

threshold_line(ax, 305)
ax.set_xlabel("PL25_M23 dose (µg/mL)", fontsize=10)
ax.set_ylabel("log$_{10}$ CFU reduction", fontsize=10)
ax.set_xticks(DOSES)
ax.set_ylim(-0.3, 5.5)
ax.legend(loc="lower right", fontsize=8.5, frameon=False, ncol=2,
          handletextpad=0.4, columnspacing=0.8, labelcolor=INK)
style(ax)
panel_letter(ax, "d")

for ext in ("pdf", "svg", "png"):
    fig6.savefig(FIGS / f"06_activity_killing.{ext}",
                 bbox_inches="tight", dpi=250, facecolor="white")
print("wrote figures/06_activity_killing.{pdf,svg,png}")


# ═════════════════════════════════════════════════════════════════════
# FIGURE 7 - biochemical properties
# ═════════════════════════════════════════════════════════════════════
fig7 = plt.figure(figsize=(13.5, 5))
gs7 = fig7.add_gridspec(1, 2, wspace=0.30)

# ── a  metal dependence ──────────────────────────────────────────────
ax = fig7.add_subplot(gs7[0, 0])
CONDITIONS = ["EDTA only", "PL25_M23", "PL25_M23+EDTA", "PL25_M23+EDTA+Zn"]
TICK_LABELS = ["EDTA only", "PL25_M23", "PL25_M23\n+ EDTA", "PL25_M23\n+ EDTA + Zn"]
x = np.arange(len(CONDITIONS))
width = 0.19

for i, strain in enumerate(STRAINS):
    sub = edta[edta.strain == strain]
    baseline = log_cfu(sub[sub.condition == "Untreated"].iloc[0]).mean()
    enzyme_alone = log_cfu(sub[sub.condition == "PL25_M23"].iloc[0])

    values = []
    for condition in CONDITIONS:
        hit = sub[sub.condition == condition]
        values.append(baseline - log_cfu(hit.iloc[0]).mean() if not hit.empty else 0)
    ax.bar(x + (i - 1.5) * width, values, width, facecolor=FILL[strain],
           edgecolor=INK, lw=0.7, hatch=HATCH[strain],
           label=LABEL[strain], zorder=3)

    # EDTA and EDTA+Zn are tested against enzyme alone, not against untreated:
    # the question is whether chelation changes killing, not whether killing occurs.
    for j, condition in enumerate(CONDITIONS[2:], start=2):
        hit = sub[sub.condition == condition]
        if hit.empty:
            continue
        p = welch(enzyme_alone, log_cfu(hit.iloc[0]))
        record("7", "a", strain, f"{condition} vs PL25_M23 alone", values[j], p)
        symbol = stars(p)
        ax.text(x[j] + (i - 1.5) * width, values[j] + 0.2, symbol,
                ha="center", va="bottom", fontsize=7 if symbol == "ns" else 9,
                fontweight="normal" if symbol == "ns" else "bold", color=INK)

ax.axhline(3, color="#888", ls=(0, (3, 3)), lw=0.7, zorder=1)
ax.set_xticks(x)
ax.set_xticklabels(TICK_LABELS, fontsize=8.5, color=INK)
ax.set_ylabel("log$_{10}$ CFU reduction", fontsize=10)
ax.set_ylim(-0.3, 5.5)
ax.legend(loc="upper left", fontsize=8.5, frameon=False, ncol=2,
          handletextpad=0.6, columnspacing=0.8, labelcolor=INK)
ax.text(0.99, 0.01, "asterisks: vs PL25_M23 alone", transform=ax.transAxes,
        fontsize=7.5, color=INK, ha="right", va="bottom", style="italic")
style(ax)
panel_letter(ax, "a")

# ── b  salt tolerance ────────────────────────────────────────────────
ax = fig7.add_subplot(gs7[0, 1])
for strain in STRAINS:
    sub = salt[salt.strain == strain]
    control = log_cfu(sub[sub.treatment == "Untreated"].iloc[0])
    baseline = None
    means, sds = [], []
    for nacl in NACL:
        treated = log_cfu(sub[(sub.treatment == "PL25_M23") & (sub.NaCl_M == nacl)].iloc[0])
        reductions = control - treated
        means.append(reductions.mean())
        sds.append(reductions.std(ddof=1))
        # Each salt condition is compared with the same strain at 0 M NaCl.
        if nacl == 0.0:
            baseline = treated
            continue
        p = welch(baseline, treated)
        record("7", "b", strain, f"{nacl} M NaCl vs 0 M", reductions.mean(), p)
        symbol = stars(p)
        ax.text(nacl, reductions.mean() + reductions.std(ddof=1) + 0.1, symbol,
                ha="center", va="bottom", fontsize=7 if symbol == "ns" else 9,
                fontweight="normal" if symbol == "ns" else "bold", color=INK)
    ax.errorbar(NACL, means, yerr=sds, marker=MARKER[strain], ms=9, capsize=3.5,
                color=INK, mfc=INK, mec=INK, mew=0.9, lw=1.6, ecolor=INK,
                label=LABEL[strain], zorder=4)

threshold_line(ax, 1.03)
ax.set_xlabel("NaCl concentration (M)", fontsize=10)
ax.set_ylabel("log$_{10}$ CFU reduction", fontsize=10)
ax.set_xticks(NACL)
ax.set_xlim(-0.08, 1.12)
ax.set_ylim(-0.3, 5.5)
ax.legend(loc="lower right", fontsize=8.5, frameon=False, ncol=2,
          handletextpad=0.4, columnspacing=0.8, labelcolor=INK)
ax.text(0.99, 0.01, "asterisks: vs 0 M NaCl", transform=ax.transAxes,
        fontsize=7.5, color=INK, ha="right", va="bottom", style="italic")
style(ax)
panel_letter(ax, "b")

for ext in ("pdf", "svg", "png"):
    fig7.savefig(FIGS / f"07_activity_biochemistry.{ext}",
                 bbox_inches="tight", dpi=250, facecolor="white")
print("wrote figures/07_activity_biochemistry.{pdf,svg,png}")

# ── Audit table of every test drawn on the two figures ───────────────
summary = pd.DataFrame(stat_rows)
summary.to_csv(DATA / "figure_stats_summary.tsv", sep="\t",
               index=False, float_format="%.4g")
print(f"wrote results/20_activity/figure_stats_summary.tsv ({len(summary)} tests)")
