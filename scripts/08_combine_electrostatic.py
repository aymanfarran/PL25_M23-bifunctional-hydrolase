"""Assemble the revised electrostatic comparison figure (Figure 8 v2).

Layout:
  2 × 5 grid of surface views  (rows = front / 180° rotation; cols = 5 proteins)
  + 1 bar chart row below (pI, net charge at pH 7, charge density %)
  + a single color key on the right of the surface grid
  + cautious caption block

All surfaces rendered with the same per-atom color rule (D/E → red, K/R → blue,
His → light pink, others → white), so panels are directly comparable.
"""
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.image import imread
import numpy as np

FIG = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"
JS  = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/16_electrostatic/pI_charge_summary.json"

with open(JS) as f: D = json.load(f)

proteins = [
    ("PL25_M23",    "PL25_M23\n(this study, AlphaFold2;\nresidues 38–370)",                  "#D7191C"),
    ("CwlP",        "CwlP-clipped\n(AlphaFold2 model; UniProt O31976;\nresidues 1373–1686)",  "#7C3AED"),
    ("CwlT",        "CwlT\n(AlphaFold2 model; UniProt P96645;\nfull-length)",                  "#FF7F0E"),
    ("Lysostaphin", "Lysostaphin\n(PDB 4QPB;\nmature)",                                        "#2C7A2C"),
    ("LytM",        "LytM\n(PDB 4ZYB;\nmature)",                                               "#E08E10"),
]

fig = plt.figure(figsize=(16, 11))
gs  = fig.add_gridspec(3, 6, height_ratios=[1.0, 1.0, 1.05],
                       width_ratios=[1, 1, 1, 1, 1, 0.22],
                       hspace=0.18, wspace=0.06)

# Row 1: front views
for i, (key, lab, col) in enumerate(proteins):
    ax = fig.add_subplot(gs[0, i])
    ax.imshow(imread(f"{FIG}/_es2_{key}_front.png"))
    ax.set_axis_off()
    ax.set_title(lab, fontsize=9.5, color=col, fontweight="bold", pad=4)
    if i == 0:
        ax.text(-0.07, 0.5, "Front view", transform=ax.transAxes,
                rotation=90, ha="right", va="center",
                fontsize=11, fontweight="bold", color="#222")
    # Highlight PL25_M23 with a red border
    if key == "PL25_M23":
        for spine in ax.spines.values():
            spine.set_visible(True); spine.set_edgecolor("#D7191C"); spine.set_linewidth(2.5)

# Row 2: 180° rotation
for i, (key, lab, col) in enumerate(proteins):
    ax = fig.add_subplot(gs[1, i])
    ax.imshow(imread(f"{FIG}/_es2_{key}_back.png"))
    ax.set_axis_off()
    if i == 0:
        ax.text(-0.07, 0.5, "180° rotation", transform=ax.transAxes,
                rotation=90, ha="right", va="center",
                fontsize=11, fontweight="bold", color="#222")
    if key == "PL25_M23":
        for spine in ax.spines.values():
            spine.set_visible(True); spine.set_edgecolor("#D7191C"); spine.set_linewidth(2.5)

# Right column: shared color key (one for both rows)
ax_key = fig.add_subplot(gs[0:2, 5])
ax_key.set_axis_off()
ax_key.text(0.5, 0.97, "Surface key", ha="center", va="top",
            fontsize=10, fontweight="bold", transform=ax_key.transAxes)
ax_key.text(0.5, 0.88,
    "Per-atom side-chain\ncharge at pH ≈ 7\n(same rule for all 5 proteins)",
    ha="center", va="top", fontsize=8, color="#444",
    transform=ax_key.transAxes, style="italic")
items = [
    ("Acidic side chain\nAsp / Glu (charge ≈ −1)", "tomato"),
    ("Basic side chain\nLys / Arg (charge ≈ +1)",  "royalblue"),
    ("Histidine\n(partial / pK ≈ 6)",              "lightpink"),
    ("Other / neutral",                            "white"),
]
y = 0.72
for label, c in items:
    ax_key.add_patch(mpatches.Rectangle((0.10, y-0.04), 0.15, 0.08,
                    facecolor=c, edgecolor="black", linewidth=0.6,
                    transform=ax_key.transAxes))
    ax_key.text(0.30, y, label, ha="left", va="center", fontsize=8,
                transform=ax_key.transAxes)
    y -= 0.13

ax_key.text(0.0, 0.05,
    "All values calculated for the modelled or experimental region used for "
    "structural comparison; PL25_M23 region = residues 38–370.",
    ha="left", va="bottom", fontsize=7.5, style="italic", color="#666",
    transform=ax_key.transAxes, wrap=True)

# Row 3: bar chart
ax3 = fig.add_subplot(gs[2, :])
labels = [k for k, _, _ in proteins]
pI    = [D["PL25_M23"]["pI"],          D["CwlP_clipped"]["pI"],          D["CwlT"]["pI"],          D["Lysostaphin"]["pI"],          D["LytM"]["pI"]]
chrg  = [D["PL25_M23"]["charge_pH7"],  D["CwlP_clipped"]["charge_pH7"],  D["CwlT"]["charge_pH7"],  D["Lysostaphin"]["charge_pH7"],  D["LytM"]["charge_pH7"]]
dens  = [D["PL25_M23"]["charge_density_pct"], D["CwlP_clipped"]["charge_density_pct"], D["CwlT"]["charge_density_pct"],
         D["Lysostaphin"]["charge_density_pct"], D["LytM"]["charge_density_pct"]]

x = np.arange(len(labels))
w = 0.27
b1 = ax3.bar(x - w,     pI,   w, color="#9CA3AF", edgecolor="black", lw=0.5, label="Predicted pI")
b2 = ax3.bar(x,         chrg, w, color="#60A5FA", edgecolor="black", lw=0.5, label="Net charge at pH 7")
b3 = ax3.bar(x + w,     dens, w, color="#F87171", edgecolor="black", lw=0.5, label="Net charge density (%)")

# Highlight PL25 column
ax3.axvspan(-0.5, 0.5, color="#FFE5E5", alpha=0.35, zorder=0)
ax3.text(0, max(pI)+1, "★ PL25_M23", ha="center", va="bottom",
         fontsize=9, fontweight="bold", color="#D7191C")

# Value labels
for bars in (b1, b2, b3):
    for r in bars:
        h = r.get_height()
        ax3.text(r.get_x()+r.get_width()/2, h + (0.5 if h>=0 else -1.2),
                 f"{h:+.1f}" if bars is not b1 else f"{h:.2f}",
                 ha="center", va="bottom" if h>=0 else "top",
                 fontsize=7.5)

ax3.axhline(0, color="#444", lw=0.6)
ax3.set_xticks(x); ax3.set_xticklabels(labels, fontsize=10)
ax3.set_ylabel("Value")
ax3.legend(loc="upper right", fontsize=8, frameon=False)
for s in ("top","right"): ax3.spines[s].set_visible(False)
ax3.set_title("Sequence-based pI, net charge, and charge density of the modelled region of each protein",
              loc="left", fontsize=10.5, fontweight="bold", pad=6)

fig.suptitle("Predicted acidic surface profile of PL25_M23 compared with reference cell-wall hydrolases",
             fontsize=12.5, fontweight="bold", y=0.995)

fig.text(0.01, 0.005,
    "Cautious interpretation: panels show a PREDICTED acidic-surface profile for PL25_M23 (pI 4.54; net charge −22.4; "
    "charge density −6.9 %) relative to four reference cell-wall hydrolases. This is consistent with — but does NOT prove — "
    "halophile-type adaptation. Confirmation requires direct biochemical tests of activity / stability under high-salt conditions. "
    "Surfaces are coloured by per-atom side-chain charge using a single rule applied identically to all five proteins; "
    "absolute vacuum-electrostatic-potential values were not used as the comparison metric because their numeric scale "
    "is protein-specific. Bar-chart values from Biopython ProteinAnalysis (Henderson–Hasselbalch, pH 7.0).",
    ha="left", va="bottom", fontsize=7.7, style="italic", color="#444", wrap=True)

fig.subplots_adjust(top=0.945, bottom=0.10, left=0.05, right=0.99, hspace=0.20, wspace=0.06)

out = f"{FIG}/08b_electrostatic_comparison_v2"
for ext in ("pdf","svg","png"):
    plt.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print(f"✓ Saved {out}.{{pdf,svg,png}}")
