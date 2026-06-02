#!/usr/bin/env python3
"""
contig_11 (74,025 bp) — putative conjugative ICE/plasmid in V. salarius PL25.
Shows all 99 Prokka CDS with HMM-confirmed mobile-element genes highlighted.
"""
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
OUT.mkdir(parents=True, exist_ok=True)

L = 74025
# Parse CDS on contig_11
cds = []
with open(GFF) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split("\t")
        if len(p) < 9 or p[0] != "contig_11" or p[2] != "CDS": continue
        m = re.search(r"locus_tag=([^;]+)", p[8])
        if not m: continue
        cds.append({"locus_tag": m.group(1),
                    "start": int(p[3]), "end": int(p[4]),
                    "strand": p[6]})
print(f"contig_11 CDS: {len(cds)}")

# Parse HMM hits — give each locus its top domain
hits = {}
with open(DOMTBL) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split()
        if len(p) < 19: continue
        target, lt, eval_, alif, alit = p[0], p[3], float(p[6]), p[17], p[18]
        if lt not in hits or eval_ < hits[lt][1]:
            hits[lt] = (target, eval_)

# Categorise CDS
CAT_COL = {
    "M23 peptidase":    "#D7191C",   # red — our enzyme of interest
    "T4SS ATPase":      "#1F77B4",   # blue
    "T4SS-DNA transfer":"#5E3C99",   # purple
    "Other CDS":        "#CCCCCC",
}
def categorise(lt):
    h = hits.get(lt, ("none", 999))
    name = h[0]
    if name == "Peptidase_M23": return "M23 peptidase"
    if name in ("AAA_10","AAA_19","AAA_23"): return "T4SS ATPase"
    if name == "T4SS-DNA_transf": return "T4SS-DNA transfer"
    return "Other CDS"

# Figure
fig, ax = plt.subplots(figsize=(16, 4))
y_fwd, y_rev = 0.45, -0.45
h = 0.30

# backbone
ax.plot([0, L], [0, 0], color="black", lw=1.0, zorder=1)

# arrows for CDS
for g in cds:
    cat = categorise(g["locus_tag"])
    col = CAT_COL[cat]
    y = y_fwd if g["strand"] == "+" else y_rev
    s, e = g["start"], g["end"]
    width = e - s
    arrow_w = min(700, width * 0.4)  # arrowhead width in bp
    body_w  = width - arrow_w
    if g["strand"] == "+":
        verts = [(s, y-h/2), (s+body_w, y-h/2), (s+body_w, y-h),
                 (e, y), (s+body_w, y+h), (s+body_w, y+h/2), (s, y+h/2), (s, y-h/2)]
    else:
        verts = [(e, y-h/2), (s+arrow_w, y-h/2), (s+arrow_w, y-h),
                 (s, y), (s+arrow_w, y+h), (s+arrow_w, y+h/2), (e, y+h/2), (e, y-h/2)]
    poly = plt.Polygon(verts, color=col, alpha=0.9 if cat != "Other CDS" else 0.5,
                       edgecolor="black" if cat != "Other CDS" else "none",
                       linewidth=0.5, zorder=3 if cat != "Other CDS" else 2)
    ax.add_patch(poly)

# Annotate the key M23 + ATPase + T4SS genes
KEY_LABELS = {
    "PL25_00076": ("Rep / Relaxase?", "#222"),
    "PL25_00078": ("M23 peptidase\n(query — our enzyme)", "#D7191C"),
    "PL25_00082": ("VirB4-like ATPase\n(T4SS coupling)", "#1F77B4"),
    "PL25_00094": ("ATPase", "#1F77B4"),
    "PL25_00148": ("T4SS-DNA transfer\n(coupling protein)", "#5E3C99"),
}
for lt, (label, col) in KEY_LABELS.items():
    g = next((x for x in cds if x["locus_tag"] == lt), None)
    if not g: continue
    mid = (g["start"] + g["end"]) / 2
    y = y_fwd if g["strand"] == "+" else y_rev
    y_label = 1.05 if y > 0 else -1.05
    ax.annotate(label,
                xy=(mid, y + (0.35 if y > 0 else -0.35)),
                xytext=(mid, y_label),
                ha="center", va="center" if y_label > 0 else "center",
                fontsize=8.5, fontweight="bold", color=col,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.8),
                zorder=10)

# Axis
ax.set_xlim(-1500, L+1500)
ax.set_ylim(-1.6, 1.6)
ax.set_xlabel("Position on contig_11 (bp)", fontsize=11)
ax.set_yticks([y_fwd, 0, y_rev])
ax.set_yticklabels(["+ strand", "", "− strand"], fontsize=10)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.grid(axis="x", alpha=0.3)

# Title & subtitle
ax.set_title("contig_11 (74,025 bp) — putative conjugative element of Virgibacillus salarius PL25\n"
             "99 CDS from independent Prokka annotation • mobile-element genes confirmed by HMM scan",
             fontsize=11, fontweight="bold")

# Legend
handles = [mpatches.Patch(color=v, label=k, alpha=0.9 if k!="Other CDS" else 0.5) for k, v in CAT_COL.items()]
ax.legend(handles=handles, loc="lower right", ncol=4, fontsize=8.5, frameon=False)

# Footnote
fig.text(0.5, 0.01,
    "HMM hits at E ≤ 1e-5 (extended ICE panel: 21 Pfam profiles). "
    "PL25_00076 (Rep/Relaxase?) is too divergent for any local Pfam HMM — IBEX BLAST vs nr will resolve identity.",
    ha="center", fontsize=8, style="italic", color="#555")

plt.tight_layout(rect=[0,0.04,1,1])
fig.savefig(OUT/"01_contig11_map.pdf", bbox_inches="tight")
fig.savefig(OUT/"01_contig11_map.svg", bbox_inches="tight")
fig.savefig(OUT/"01_contig11_map.png", dpi=200, bbox_inches="tight")
print(f"✓ {OUT/'01_contig11_map.pdf|svg|png'}")
