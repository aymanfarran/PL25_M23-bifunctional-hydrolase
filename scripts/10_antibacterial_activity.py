#!/usr/bin/env python3
"""Wet-lab antibacterial activity of LysPep_M23 (= PL25_M23).
Computes mean ± SD, log10 reduction, % killing, Welch's t-test (on log10 CFU),
and Mann-Whitney U test. Renders a publication-quality bar chart.

Saves:
  figures/09_antibacterial_activity.{pdf,svg,png}
  results/20_activity/activity_summary.tsv
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
OUT  = ROOT/"figures"
RES  = ROOT/"results/20_activity"; RES.mkdir(parents=True, exist_ok=True)

# ── Raw data (CFU/mL) — updated 2026-06-02 ────────────────
# Assay: PL25_M23 (100 µg/mL) vs OD600-normalised bacteria, 2 h at 37 °C in PBS,
#        serial dilution + spotting on LB agar; 3 technical replicates per condition.
data = {
    "Bacillus subtilis":              {"Neg":[5e8, 4e8, 5e8],    "M23":[4e5, 7e5, 2e5]},
    "Staphylococcus warneri":         {"Neg":[5e8, 3e8, 5e8],    "M23":[2e6, 8e5, 4e6]},
    "MRSA":                           {"Neg":[6e8, 3e8, 5e8],    "M23":[4e6, 8e5, 2e6]},
    "Virgibacillus salarius PL25":    {"Neg":[4e8, 5e8, 3e8],    "M23":[7e4, 3e4, 7e3]},
}

# ── Stats ────────────────────────────────────────────────
rows = []
for strain, d in data.items():
    ctrl     = np.array(d["Neg"])
    treated  = np.array(d["M23"])
    log_c    = np.log10(ctrl)
    log_t    = np.log10(treated)
    mean_c   = ctrl.mean();   sd_c = ctrl.std(ddof=1)
    mean_t   = treated.mean();sd_t = treated.std(ddof=1)
    log_red  = log_c.mean() - log_t.mean()
    pct_kill = (1 - mean_t/mean_c) * 100

    # Welch's t-test on log10-transformed CFU (standard for log-normal microbial data)
    t_stat, p_t = stats.ttest_ind(log_c, log_t, equal_var=False)
    # Mann-Whitney U (non-parametric back-up; with n=3 min p ≈ 0.05)
    u_stat, p_u = stats.mannwhitneyu(ctrl, treated, alternative="two-sided")

    rows.append({
        "Strain": strain,
        "Ctrl_mean_CFU_per_mL": mean_c,
        "Ctrl_SD": sd_c,
        "M23_mean_CFU_per_mL": mean_t,
        "M23_SD": sd_t,
        "Log10_reduction": log_red,
        "Percent_killing": pct_kill,
        "Welch_t":   t_stat,
        "p_Welch":   p_t,
        "MW_U":      u_stat,
        "p_MW":      p_u,
        "n": 3,
    })

df = pd.DataFrame(rows)
print(df.to_string(index=False))
df.to_csv(RES/"activity_summary.tsv", sep="\t", index=False, float_format="%.4g")
print(f"\n✓ Wrote {RES}/activity_summary.tsv")

# ── Figure: bar chart of log10 CFU per strain × treatment + replicate dots + p-values ────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                               gridspec_kw={"width_ratios":[1.4, 1.0]})

strains = list(data.keys())
short_lab = {
    "Bacillus subtilis":              "$\\it{B.\\ subtilis}$",
    "Staphylococcus warneri":         "$\\it{S.\\ warneri}$",
    "MRSA":                           "MRSA",
    "Virgibacillus salarius PL25":    "$\\it{V.\\ salarius}$\nPL25",
}
x = np.arange(len(strains))
w = 0.38

# --- Panel A: log10 CFU/mL bars with dots and p-values ---
log_ctrl_means = [np.mean(np.log10(data[s]["Neg"])) for s in strains]
log_treat_means= [np.mean(np.log10(data[s]["M23"])) for s in strains]
log_ctrl_sds   = [np.std(np.log10(data[s]["Neg"]), ddof=1) for s in strains]
log_treat_sds  = [np.std(np.log10(data[s]["M23"]), ddof=1) for s in strains]

b1 = ax1.bar(x-w/2, log_ctrl_means,  w, yerr=log_ctrl_sds,  capsize=4,
             color="#9CA3AF", edgecolor="black", lw=0.7, label="Negative control")
b2 = ax1.bar(x+w/2, log_treat_means, w, yerr=log_treat_sds, capsize=4,
             color="#D7191C", edgecolor="black", lw=0.7, label="PL25_M23 (100 µg/mL)")

# Replicate dots overlay
for i, s in enumerate(strains):
    for v in data[s]["Neg"]:
        ax1.plot(i-w/2, np.log10(v), "o", color="black", markersize=4, alpha=0.55, mfc="white")
    for v in data[s]["M23"]:
        ax1.plot(i+w/2, np.log10(v), "o", color="black", markersize=4, alpha=0.55, mfc="white")

# Significance annotations
for i, s in enumerate(strains):
    p = df.loc[df.Strain==s, "p_Welch"].iloc[0]
    if   p < 0.001: tag = "***"
    elif p < 0.01:  tag = "**"
    elif p < 0.05:  tag = "*"
    else:           tag = "ns"
    h = max(log_ctrl_means[i], log_treat_means[i]) + max(log_ctrl_sds[i], log_treat_sds[i]) + 0.3
    ax1.plot([i-w/2, i-w/2, i+w/2, i+w/2], [h-0.1, h, h, h-0.1], "k-", lw=0.7)
    ax1.text(i, h+0.05, f"{tag}\np = {p:.2g}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

ax1.set_xticks(x)
ax1.set_xticklabels([short_lab[s] for s in strains], fontsize=10)
ax1.set_ylabel("log₁₀ CFU mL⁻¹", fontsize=11)
ax1.set_ylim(2, 11.5)
ax1.set_title("a   Bacterial viability after 2 h exposure to PL25_M23 (100 µg/mL, 37 °C, PBS)",
              loc="left", fontsize=11, fontweight="bold")
ax1.legend(loc="upper right", frameon=False, fontsize=9)
for sp in ("top","right"): ax1.spines[sp].set_visible(False)

# --- Panel B: log10 reduction bar with %-killing labels ---
log_reds = [df.loc[df.Strain==s, "Log10_reduction"].iloc[0] for s in strains]
pct_kls  = [df.loc[df.Strain==s, "Percent_killing"].iloc[0] for s in strains]
cols     = ["#E5E7EB","#9CA3AF","#D97706","#D7191C"]   # neutral palette; emphasise MRSA + self

bars = ax2.bar(x, log_reds, 0.6, color=cols, edgecolor="black", lw=0.7)
for r, pct, l in zip(bars, pct_kls, log_reds):
    ax2.text(r.get_x()+r.get_width()/2, l+0.08,
             f"{l:.2f} log\n({pct:.2f} % killed)",
             ha="center", va="bottom", fontsize=9, fontweight="bold")

ax2.set_xticks(x)
ax2.set_xticklabels([short_lab[s] for s in strains], fontsize=10)
ax2.set_ylabel("log₁₀ reduction (vs. negative control)", fontsize=11)
ax2.set_ylim(0, max(log_reds)+1.2)
ax2.set_title("b   Magnitude of bactericidal effect", loc="left", fontsize=11, fontweight="bold")
for sp in ("top","right"): ax2.spines[sp].set_visible(False)
ax2.axhline(2, ls="--", lw=0.8, color="#888")
ax2.text(len(strains)-0.5, 2.05, "2-log = 99 % killing", ha="right", va="bottom",
         fontsize=8, style="italic", color="#666")

fig.suptitle("PL25_M23 shows bactericidal activity against Gram-positive bacteria, including MRSA",
             fontsize=12.5, fontweight="bold", y=1.005)

fig.text(0.01, 0.005,
    "Assay: OD₆₀₀-normalised bacterial suspensions were incubated with PL25_M23 at 100 µg/mL or with PBS "
    "(negative control) for 2 h at 37 °C, serially diluted, and spotted on LB agar to enumerate surviving CFU. "
    "Bars in (a) show mean log₁₀ CFU mL⁻¹ ± SD with individual replicates as black dots (n = 3 per condition). "
    "Significance: Welch's t-test on log₁₀-transformed CFU/mL data; *p < 0.05; **p < 0.01; ***p < 0.001; ns = not significant. "
    "Bars in (b) show the mean log₁₀ reduction in viable counts and the corresponding percentage killing. "
    "Data are from a single experiment with three technical replicates; biological-replicate confirmation is recommended.",
    ha="left", va="bottom", fontsize=8, style="italic", color="#555")

plt.tight_layout(rect=[0, 0.05, 1, 0.97])
out = OUT/"09_antibacterial_activity"
for ext in ("pdf","svg","png"):
    fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print(f"✓ Saved {out}.{{pdf,svg,png}}")
