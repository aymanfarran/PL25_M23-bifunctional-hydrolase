"""Assemble the 4-panel 3D structure figure (v2) with full-length numbering."""
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.patches import Patch

FIG = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"

panels = [
    ("_panel_A_front_v2.png",      "A — Front view",
     "Slt35-like glycan hydrolase (blue, N-terminal) + M23 peptidase (red, C-terminal)"),
    ("_panel_B_side_v2.png",       "B — 90° rotation",
     "Two-domain architecture clearly separated"),
    ("_panel_C_back_v2.png",       "C — 180° back view",
     "Linker (grey) bridges the two catalytic domains"),
    ("_panel_D_activesite_v2.png", "D — M23 putative Zn²⁺-binding site",
     "H264, D268, H346, H348  (HxxxD + HxH motifs, full-length numbering)"),
]

fig, axes = plt.subplots(2, 2, figsize=(13, 11))
for ax, (img, title, sub) in zip(axes.ravel(), panels):
    ax.imshow(imread(f"{FIG}/{img}"))
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=6)
    ax.text(0.5, -0.05, sub, transform=ax.transAxes, ha="center", va="top",
            fontsize=9, color="#333", style="italic")
    ax.set_axis_off()

fig.suptitle("PL25_M23 — bifunctional cell-wall hydrolase, CwlP-like architecture (AlphaFold2)\n"
             "Full-length 370 aa  •  modelled region (residues 38–370)  •  mean pLDDT = 88.7",
             fontsize=12, fontweight="bold", y=0.985)

legend = [
    Patch(color="#2C7A2C", label="Slt35-like glycan hydrolase fold (residues 70–222)"),
    Patch(color="#E9E9A6", label="Linker region (residues 223–262)"),
    Patch(color="#B22222", label="M23 peptidase (residues 263–360; PF01551)"),
    Patch(color="#FFA500", label="Putative Zn²⁺-binding residues (H264, D268, H346, H348)"),
]
fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.02),
           ncol=2, fontsize=9, frameon=True, framealpha=0.95)

fig.subplots_adjust(top=0.92, bottom=0.07, left=0.02, right=0.98, wspace=0.05, hspace=0.18)

out = f"{FIG}/04c_M23_bifunctional_panel_v2"
for ext in ("pdf","svg","png"):
    plt.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print(f"✓ Saved {out}.{{pdf,svg,png}}")
