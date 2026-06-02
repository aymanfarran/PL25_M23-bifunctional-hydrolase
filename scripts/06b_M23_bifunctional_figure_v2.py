#!/usr/bin/env python3
"""
PL25 M23 — BIFUNCTIONAL architecture figure (v2).

CHANGE vs 03b:
- Removed any implication of a Pfam PF01464 (SLT) assignment for the
  N-terminal domain. InterProScan returns NO hit for PF01464 on PL25_M23;
  the lysozyme-/Slt35-like fold is detected only at the FOLD level by
  CDD cd13399, SUPERFAMILY SSF53955 and Gene3D G3DSA:1.10.530.10.
- Label now reads: "Slt35-/Lysozyme-like fold (CDD cd13399; SSF53955;
  Gene3D)  —  fold-level assignment; below Pfam PF01464 HMM threshold".
- Title updated: "CwlT-like" → "CwlP-like" (CwlT is not in the M23 family).

Saves to figures/03c_M23_bifunctional_architecture_v2.{pdf,svg,png}.
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

PROJ = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
PDB  = PROJ/"results/08_colabfold/PL25_M23_mature/PL25_M23_mature_PL25_M23_mature_form__signal_peptide_removed__residues_38-370__unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"

plddt = {}
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
seq = {}
with open(PDB) as f:
    for line in f:
        if not line.startswith("ATOM"): continue
        if line[12:16].strip() != "CA": continue
        r = int(line[22:26])
        plddt[r] = float(line[60:66])
        seq[r] = AA3.get(line[17:20].strip(), "X")

DOM_SLT35  = (33, 185)
DOM_LINKER = (186, 225)
DOM_M23    = (226, 323)
ZN_SITE    = [227, 231, 309, 311]

sequence = "".join(seq[r] for r in sorted(seq.keys()))
KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"E":-3.5,"Q":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
w = 19
hydro = []
for i in range(len(sequence)-w+1):
    hydro.append(sum(KD.get(a,0) for a in sequence[i:i+w])/w)

fig = plt.figure(figsize=(15, 6.8))
gs = fig.add_gridspec(3, 1, height_ratios=[2.2, 0.6, 0.95], hspace=0.20)

# pLDDT
ax = fig.add_subplot(gs[0])
xs = list(plddt.keys()); ys = [plddt[r] for r in xs]
colors = ["#0053D6" if v>=90 else "#65CBF3" if v>=70 else "#FFDB13" if v>=50 else "#FF7D45" for v in ys]
ax.bar(xs, ys, width=1.0, color=colors, edgecolor="none")
ax.axhline(90, ls="--", color="#0053D6", lw=0.5, alpha=0.5)
ax.axhline(70, ls="--", color="#65CBF3", lw=0.5, alpha=0.5)
for r in ZN_SITE:
    if r in plddt:
        ax.annotate(f"{seq[r]}{r}", xy=(r, plddt[r]), xytext=(r, 103),
                    ha="center", va="bottom", fontsize=7.5, color="#E65100", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#E65100", lw=0.8))
ax.set_ylabel("pLDDT", fontsize=10)
ax.set_ylim(0, 100); ax.set_xlim(0, max(xs)+5); ax.set_xticklabels([])
mean_plddt = np.mean(list(plddt.values()))
ax.text(0.99, 0.96, f"Mean pLDDT = {mean_plddt:.1f}", transform=ax.transAxes, ha="right", va="top",
        fontsize=10, fontweight="bold", bbox=dict(facecolor="white", edgecolor="#888", alpha=0.9))
ax.set_title("PL25 M23 — bifunctional cell-wall hydrolase (Slt35-/Lysozyme-like fold + Pfam M23, CwlP-like architecture)",
             fontsize=12, fontweight="bold")

# Hydropathy
ax2 = fig.add_subplot(gs[1], sharex=ax)
ax2.plot(range(w//2+1, w//2+1+len(hydro)), hydro, color="#444", lw=1)
ax2.fill_between(range(w//2+1, w//2+1+len(hydro)), hydro, 0,
                 where=[h>0 for h in hydro], color="#FFA726", alpha=0.5)
ax2.axhline(0, color="#888", lw=0.5)
ax2.axhline(1.6, color="red", lw=0.6, ls="--", alpha=0.5)
ax2.set_ylabel("KD", fontsize=9); ax2.set_xticklabels([])
for s in ("top","right"): ax2.spines[s].set_visible(False)

# Domain bar — CORRECTED LABEL
ax3 = fig.add_subplot(gs[2], sharex=ax)
domain_h = 0.45
ax3.add_patch(mpatches.Rectangle((DOM_SLT35[0], 0), DOM_SLT35[1]-DOM_SLT35[0], domain_h,
              facecolor="#2196F3", alpha=0.92, edgecolor="black", linewidth=0.5))
ax3.text((DOM_SLT35[0]+DOM_SLT35[1])/2, domain_h/2,
         "Slt35-/Lysozyme-like fold\n(CDD cd13399; SSF53955; Gene3D)\nGLYCAN HYDROLASE — fold-level assignment",
         ha="center", va="center", fontsize=8, color="white", fontweight="bold")
ax3.add_patch(mpatches.Rectangle((DOM_LINKER[0], 0), DOM_LINKER[1]-DOM_LINKER[0], domain_h,
              facecolor="#E0E0E0", alpha=0.7, edgecolor="black", linewidth=0.5))
ax3.text((DOM_LINKER[0]+DOM_LINKER[1])/2, domain_h/2,
         "linker", ha="center", va="center", fontsize=7, color="black", style="italic")
ax3.add_patch(mpatches.Rectangle((DOM_M23[0], 0), DOM_M23[1]-DOM_M23[0], domain_h,
              facecolor="#D7191C", alpha=0.92, edgecolor="black", linewidth=0.5))
ax3.text((DOM_M23[0]+DOM_M23[1])/2, domain_h/2,
         "M23 peptidase\n(Pfam PF01551, Zn²⁺)\nPEPTIDE HYDROLASE",
         ha="center", va="center", fontsize=8, color="white", fontweight="bold")
for r in ZN_SITE:
    ax3.plot([r,r], [-0.08, domain_h+0.08], color="orange", lw=2, zorder=10)
    ax3.scatter([r], [-0.12], marker="v", s=40, color="orange",
                edgecolor="black", linewidths=0.5, zorder=11)
ax3.text(np.mean(ZN_SITE), -0.22, "Zn²⁺ active-site residues  (canonical HxxxD + HxH)",
         ha="center", va="top", fontsize=8, color="#E65100", fontweight="bold")
ax3.set_xlim(0, max(xs)+5)
ax3.set_ylim(-0.35, domain_h+0.10)
ax3.set_yticks([])
ax3.set_xlabel("Residue (mature form numbering, ∆ signal anchor)")
for s in ("top","right","left"): ax3.spines[s].set_visible(False)

fig.text(0.01, 0.012,
    "N-terminal TM anchor (residues 19–46 of full-length protein) removed from mature form. "
    "SignalP6 predicts OTHER = 1.000 (no cleavable signal peptide). "
    "N-terminal domain is detected at the FOLD level by CDD cd13399 (E = 1.1×10⁻²⁷), "
    "SUPERFAMILY SSF53955 (E = 1.5×10⁻²⁵) and Gene3D, but is BELOW the Pfam PF01464 HMM detection threshold. "
    "The CwlP-like architectural assignment therefore rests on fold homology rather than Pfam-family membership.",
    fontsize=8, style="italic", color="#444", ha="left", wrap=True)

items = [
    mpatches.Patch(color="#2196F3", label="Slt35-/Lysozyme-like glycoside hydrolase (fold-level)"),
    mpatches.Patch(color="#D7191C", label="M23 metallopeptidase (Pfam PF01551)"),
    mpatches.Patch(color="#FFA726", label="N-terminal hydrophobic TM anchor (full-length only)"),
    mpatches.Patch(color="orange",  label="Zn²⁺ catalytic residues"),
]
fig.legend(handles=items, loc="lower right", bbox_to_anchor=(0.99, 0.005),
           ncol=2, fontsize=8.5, frameon=True, framealpha=0.95)

plt.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(PROJ/"figures/03c_M23_bifunctional_architecture_v2.pdf", bbox_inches="tight")
fig.savefig(PROJ/"figures/03c_M23_bifunctional_architecture_v2.svg", bbox_inches="tight")
fig.savefig(PROJ/"figures/03c_M23_bifunctional_architecture_v2.png", dpi=200, bbox_inches="tight")
print("✓ figures/03c_M23_bifunctional_architecture_v2.{pdf,svg,png}")
