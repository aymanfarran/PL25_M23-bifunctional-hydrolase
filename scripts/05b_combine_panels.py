"""Combine the 5 per-protein active-site close-ups + overlay into a single
publication-style figure."""
import matplotlib.pyplot as plt
from matplotlib.image import imread

FIG = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"
panels = [
    ("PL25", "PL25_M23\n(this study, AF2)", "#D7191C"),
    ("LYSO", "Lysostaphin\n(PDB 4QPB)",      "#2C7A2C"),
    ("LYTM", "LytM\n(PDB 4ZYB)",             "#E08E10"),
    ("ALE1", "ALE-1\n(AF2, O05156)",         "#1F77B4"),
    ("CWLP", "CwlP-clipped\n(AF2)",          "#7C3AED"),
]

fig = plt.figure(figsize=(16, 7.5))
gs  = fig.add_gridspec(2, 5, height_ratios=[1.6, 1.0], hspace=0.18, wspace=0.06)

# Top row: large overlay spanning all 5 columns
ax_top = fig.add_subplot(gs[0, :])
ax_top.imshow(imread(f"{FIG}/09_active_site_overlay.png"))
ax_top.set_axis_off()
ax_top.set_title("a   Overlay of M23 catalytic-site residues — all five proteins superposed on PL25_M23 (residues 205–330)",
                 loc="left", fontsize=12, fontweight="bold", pad=8)

# Bottom row: 5 per-protein close-ups
for i, (k, lab, col) in enumerate(panels):
    ax = fig.add_subplot(gs[1, i])
    ax.imshow(imread(f"{FIG}/09_active_site_{k}.png"))
    ax.set_axis_off()
    ax.set_title(lab, fontsize=10, color=col, fontweight="bold", pad=4)

fig.text(0.5, 0.015,
    "b   Per-protein close-ups of the canonical HxxxD + HxH Zn²⁺-binding residues. "
    "Side chains shown as sticks; residue numbering as in the respective parent structure.",
    ha="center", fontsize=9.5, style="italic", color="#444")

out = f"{FIG}/09_active_site_comparison"
for ext in ("pdf","svg","png"):
    plt.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print(f"✓ Saved {out}.{{pdf,svg,png}}")
