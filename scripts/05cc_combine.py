"""Assemble the 3-panel Figure 5 (revised) using the corrected CwlT structure
(AF2 of UniProt P96645). Saves figures/05cc_M23_3way_alignment.{pdf,svg,png}."""
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.patches import Patch

FIG = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"

panels = [
    ("_05cc_panelA_M23only.png",
     "a   C-terminal peptidase domain\n(domain-restricted alignment)",
     "PL25 vs CwlP M23 : 0.69 Å (541 atoms)\nPL25 vs CwlT NlpC/P60 : 3.23 Å (81 atoms — different family)"),
    ("_05cc_panelB_Nterm.png",
     "b   N-terminal glycan-hydrolase domain\n(domain-restricted alignment)",
     "PL25 vs CwlP SLT : 3.31 Å (542 atoms)\nPL25 vs CwlT GH25 : 10.43 Å (579 atoms — distinct fold)"),
    ("_05cc_panelC_full.png",
     "c   Full predicted two-domain region\n(unrestricted alignment)",
     "PL25 vs CwlP : 16.75 Å (1362 atoms)\nPL25 vs CwlT : 10.43 Å (579 atoms)"),
]

fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.2))
for ax, (img, title, sub) in zip(axes, panels):
    ax.imshow(imread(f"{FIG}/{img}"))
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", pad=6)
    ax.text(0.5, -0.08, sub, transform=ax.transAxes, ha="center", va="top",
            fontsize=9, color="#333", style="italic")
    ax.set_axis_off()

# Legend
handles = [
    Patch(color="#D7191C", label="PL25_M23 (this study, AF2) — predicted SLT + M23"),
    Patch(color="#7C3AED", label="CwlP-clipped (AF2 of O31976) — SLT + M23"),
    Patch(color="#FF7F0E", label="CwlT (AF2 of P96645) — GH25 + NlpC/P60   [different family — shown for contrast]"),
]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05),
           ncol=3, fontsize=9, frameon=True, framealpha=0.9)

fig.subplots_adjust(top=0.92, bottom=0.18, left=0.02, right=0.98, wspace=0.05)

out = f"{FIG}/05cc_M23_3way_alignment"
for ext in ("pdf","svg","png"):
    plt.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print(f"✓ Saved {out}.{{pdf,svg,png}}")
