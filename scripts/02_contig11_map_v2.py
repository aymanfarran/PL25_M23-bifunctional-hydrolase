#!/usr/bin/env python3
"""contig_11 (74,025 bp) — V. salarius PL25 conjugative plasmid map (v2 — with BLAST annotations)."""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

PROJ = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
GFF  = PROJ/"results/01_annotation/prokka/PL25_indep.gff"
DOMTBL = PROJ/"results/03_M23_ICE_locus/contig_11_extended.domtbl"
OUT  = PROJ/"figures"

L = 74025

# Parse all 99 CDS on contig_11
cds = []
with open(GFF) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split("\t")
        if len(p) < 9 or p[0] != "contig_11" or p[2] != "CDS": continue
        m = re.search(r"locus_tag=([^;]+)", p[8])
        if not m: continue
        cds.append({"locus_tag": m.group(1), "start": int(p[3]), "end": int(p[4]), "strand": p[6]})

# HMM hits
hits = {}
with open(DOMTBL) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split()
        if len(p) < 19: continue
        target, lt, eval_ = p[0], p[3], float(p[6])
        if lt not in hits or eval_ < hits[lt][1]:
            hits[lt] = (target, eval_)

# Known functional annotations from BLAST + HMM (manually curated from results)
ANNOT = {
    "PL25_00076": ("Replication-Relaxation", "#9C27B0"),
    "PL25_00077": ("RusA resolvase",         "#673AB7"),
    "PL25_00078": ("M23 peptidase (lysin)",  "#D7191C"),
    "PL25_00082": ("VirB4-like ATPase",      "#1F77B4"),
    "PL25_00094": ("ATPase (helicase)",      "#1F77B4"),
    "PL25_00148": ("T4SS-DNA transfer",      "#5E3C99"),
}

def cat_of(lt):
    if lt in ANNOT: return "function"
    h = hits.get(lt, ("",999))[0]
    if h in ("Peptidase_M23","AAA_10","AAA_19","AAA_23","T4SS-DNA_transf","Phage_integrase"): return "mge"
    return "hypothetical"

# Figure
fig, ax = plt.subplots(figsize=(16, 4.2))
y_fwd, y_rev = 0.45, -0.45
h = 0.30

ax.plot([0, L], [0, 0], color="black", lw=1.0, zorder=1)

for g in cds:
    cat = cat_of(g["locus_tag"])
    col = ANNOT[g["locus_tag"]][1] if cat == "function" else ("#FFA726" if cat == "mge" else "#CCCCCC")
    y = y_fwd if g["strand"] == "+" else y_rev
    s, e = g["start"], g["end"]
    width = e - s
    arrow_w = min(700, width * 0.4)
    body_w  = width - arrow_w
    if g["strand"] == "+":
        verts = [(s, y-h/2), (s+body_w, y-h/2), (s+body_w, y-h),
                 (e, y), (s+body_w, y+h), (s+body_w, y+h/2), (s, y+h/2)]
    else:
        verts = [(e, y-h/2), (s+arrow_w, y-h/2), (s+arrow_w, y-h),
                 (s, y), (s+arrow_w, y+h), (s+arrow_w, y+h/2), (e, y+h/2)]
    poly = plt.Polygon(verts, color=col, alpha=0.92 if cat != "hypothetical" else 0.45,
                       edgecolor="black" if cat == "function" else ("#666" if cat == "mge" else "none"),
                       linewidth=0.6 if cat != "hypothetical" else 0, zorder=4 if cat == "function" else (3 if cat == "mge" else 2))
    ax.add_patch(poly)

# Labels for annotated genes (curated)
LABEL_OFFSETS = {  # (vertical offset multiplier — 1.0 = standard above)
    "PL25_00076": (0.90, "left"),
    "PL25_00077": (1.45, "center"),
    "PL25_00078": (0.90, "right"),
    "PL25_00082": (0.90, "center"),
    "PL25_00094": (0.90, "center"),
    "PL25_00148": (0.90, "center"),
}
for lt, (label, col) in ANNOT.items():
    g = next((x for x in cds if x["locus_tag"] == lt), None)
    if not g: continue
    mid = (g["start"] + g["end"]) / 2
    y = y_fwd if g["strand"] == "+" else y_rev
    yoff, ha = LABEL_OFFSETS.get(lt, (0.90, "center"))
    y_label = (1.15 + yoff*0.0) if y > 0 else -(1.10 + yoff*0.05)
    ax.annotate(f"{lt}\n{label}",
                xy=(mid, y + (0.35 if y > 0 else -0.35)),
                xytext=(mid, y_label),
                ha=ha, va="bottom" if y_label > 0 else "top",
                fontsize=8.5, fontweight="bold", color=col,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.8),
                zorder=10)

# Axis
ax.set_xlim(-1500, L+1500)
ax.set_ylim(-2.0, 2.0)
ax.set_xlabel("Position on contig_11 (bp)", fontsize=11)
ax.set_yticks([y_fwd, 0, y_rev])
ax.set_yticklabels(["+ strand", "", "− strand"], fontsize=10)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.grid(axis="x", alpha=0.3)

ax.set_title("contig_11 (74,025 bp) — Virgibacillus salarius PL25 conjugative plasmid\n"
             "geNomad plasmid score = 0.989  •  99 CDS  •  closest characterised relative: Ornithinibacillus contaminans (57–67% id)",
             fontsize=11.5, fontweight="bold")

# Legend
legend_items = [
    mpatches.Patch(color="#D7191C", label="M23 peptidase (this paper)"),
    mpatches.Patch(color="#9C27B0", label="Replication / Relaxation"),
    mpatches.Patch(color="#673AB7", label="RusA Holliday resolvase"),
    mpatches.Patch(color="#1F77B4", label="T4SS ATPase"),
    mpatches.Patch(color="#5E3C99", label="T4SS DNA-transfer"),
    mpatches.Patch(color="#FFA726", label="Other MGE gene (HMM)", alpha=0.92),
    mpatches.Patch(color="#CCCCCC", label="hypothetical CDS", alpha=0.45),
]
ax.legend(handles=legend_items, loc="lower right", ncol=4, fontsize=8, frameon=False)

# Footnote
fig.text(0.5, 0.01,
    "Independent Prokka annotation • HMM panel: 21 Pfam profiles (T4SS / relaxase / integrase / M23) • "
    "BLAST vs NCBI nr confirms conjugation machinery on minus strand.",
    ha="center", fontsize=8, style="italic", color="#555")

plt.tight_layout(rect=[0,0.04,1,1])
fig.savefig(OUT/"01_contig11_map_v2.pdf", bbox_inches="tight")
fig.savefig(OUT/"01_contig11_map_v2.svg", bbox_inches="tight")
fig.savefig(OUT/"01_contig11_map_v2.png", dpi=200, bbox_inches="tight")
print(f"✓ Saved 01_contig11_map_v2.{{pdf,svg,png}}")
