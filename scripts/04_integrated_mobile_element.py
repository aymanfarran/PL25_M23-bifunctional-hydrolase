#!/usr/bin/env python3
"""
Integrated 'mobile element' figure — combines contig_11 plasmid features in one panel:
  • Gene map (contig_11) with mobility genes colour-coded
  • M23 location starred
  • GC deviation profile (sliding window vs chromosome baseline)
  • Synteny ribbons to closest characterised relative (O. contaminans)
  • Gene map (O. contaminans contig NZ_LDPV02000047.1)
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from Bio import SeqIO
import numpy as np

ROOT = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
SYNT = ROOT/"results/14_synteny"

# Parse GBK files
def parse_gbk(path):
    rec = next(SeqIO.parse(str(path),"genbank"))
    genes = []
    for f in rec.features:
        if f.type != "CDS": continue
        q = f.qualifiers
        lt = q.get("protein_id", q.get("locus_tag",["?"]))[0]
        prod = q.get("product",["hypothetical"])[0]
        try:
            s = int(f.location.start)
            e = int(f.location.end)
        except: continue
        strand = "+" if f.location.strand == 1 else "-"
        genes.append({"name": lt, "start": s, "end": e, "strand": strand, "product": prod})
    try: seq = str(rec.seq)
    except: seq = ""
    return genes, len(rec.seq), seq

pl25_genes, pl25_L, pl25_seq  = parse_gbk(SYNT/"PL25_contig11.gbk")
oc_genes,   oc_L,   _   = parse_gbk(SYNT/"NZ_LDPV02000047.1.gbk")

# PL25 chromosome mean GC (computed in WP2 analysis)
CHROM_GC = 36.81

# Compute GC deviation for contig_11
def gc_window(seq, window=1500, step=500):
    out = []
    for i in range(0, len(seq)-window, step):
        chunk = seq[i:i+window]
        gc = (chunk.count("G")+chunk.count("C"))/len(chunk)*100
        out.append((i+window/2, gc))
    return np.array(out)

gc_pl25 = gc_window(pl25_seq)

# Key mobility genes (from BLAST + HMM)
KEY_GENES = {
    "PL25_00076": ("Replication-Relaxation",         "#9C27B0"),
    "PL25_00077": ("RusA resolvase",                 "#673AB7"),
    "PL25_00078": ("M23 peptidase (this study)",     "#D7191C"),
    "PL25_00082": ("VirB4-like ATPase",              "#1976D2"),
    "PL25_00094": ("RecD2-like helicase",            "#42A5F5"),
    "PL25_00148": ("T4SS-DNA transfer",              "#5E3C99"),
}

# Read clinker alignments (PL25 vs O. contaminans only)
alignments = []
in_section = False; cur_q, cur_t = None, None
with open(SYNT/"PL25_vs_relatives_alignments.csv") as f:
    for raw in f:
        line = raw.rstrip()
        if not line.strip():
            in_section = False; continue
        if line.startswith("---"): continue
        if " vs " in line and not line.startswith("Query"):
            parts = line.split(" vs ")
            cur_q, cur_t = parts[0].strip(), parts[1].strip()
            in_section = False; continue
        parts = line.split()
        if len(parts) < 4: continue
        if parts[0]=="Query" and parts[1]=="Target":
            in_section=True; continue
        if in_section and cur_q=="PL25_contig11" and cur_t=="NZ_LDPV02000047.1":
            try:
                q, t, ident = parts[0], parts[1], float(parts[2])
                alignments.append((q, t, ident))
            except: pass
print(f"PL25 ↔ O. contaminans alignments: {len(alignments)}")

# Gene lookup
pl25_pos = {g["name"]: g for g in pl25_genes}
oc_pos   = {g["name"]: g for g in oc_genes}

# ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 9))
gs = GridSpec(5, 1, height_ratios=[1.5, 0.9, 1.8, 0.9, 1.5], hspace=0.10)

# ── Panel A: PL25 contig_11 gene map ───────────────────────────
ax1 = fig.add_subplot(gs[0])

def draw_gene(ax, s, e, y, strand, col, height=0.45, edgecolor=None, alpha=1.0):
    arrow_w = min(800, (e-s)*0.35)
    body_w = (e-s) - arrow_w
    if strand=="+":
        xs = [s, s+body_w, s+body_w, e, s+body_w, s+body_w, s]
        ys = [y-height/2, y-height/2, y-height, y, y+height, y+height/2, y+height/2]
    else:
        xs = [e, s+arrow_w, s+arrow_w, s, s+arrow_w, s+arrow_w, e]
        ys = [y-height/2, y-height/2, y-height, y, y+height, y+height/2, y+height/2]
    ax.fill(xs, ys, color=col, alpha=alpha, lw=0.4,
            edgecolor=edgecolor if edgecolor else "black", zorder=4)

# Plus strand at y=0.5, minus strand at y=-0.5
ax1.plot([0, pl25_L], [0, 0], color="#444", lw=1.0)
for g in pl25_genes:
    y = 0.5 if g["strand"]=="+" else -0.5
    col = KEY_GENES.get(g["name"], (None, "#BBB"))[1]
    alpha = 1.0 if g["name"] in KEY_GENES else 0.6
    edge = "black" if g["name"] in KEY_GENES else "none"
    draw_gene(ax1, g["start"], g["end"], y, g["strand"], col, height=0.40,
              edgecolor=edge, alpha=alpha)
# Star + label on M23
m23 = pl25_pos.get("PL25_00078")
if m23:
    mid = (m23["start"]+m23["end"])/2
    ax1.scatter([mid], [0.5 if m23["strand"]=="+" else -0.5], marker="*",
                s=300, color="#FFD700", edgecolor="black", linewidths=1, zorder=10)
    ax1.annotate(f"M23 peptidase\n({m23['name']})",
                 xy=(mid, 0.5 if m23["strand"]=="+" else -0.5),
                 xytext=(mid, 1.6),
                 ha="center", va="bottom", fontsize=10, fontweight="bold", color="#D7191C",
                 arrowprops=dict(arrowstyle="-|>", color="#D7191C", lw=1.2), zorder=11)
# Labels for other key genes — stagger to avoid overlap
LABEL_Y = {
    "PL25_00076": -1.7,
    "PL25_00077": -2.0,
    "PL25_00082": -1.7,
    "PL25_00094": -2.0,
    "PL25_00148": -1.7,
}
LABEL_NUDGE = {  # x-offset in bp to spread out labels
    "PL25_00076": -2500,
    "PL25_00077": +2500,
    "PL25_00082":  0,
    "PL25_00094":  0,
    "PL25_00148":  0,
}
for lt, (label, col) in KEY_GENES.items():
    if lt == "PL25_00078": continue
    g = pl25_pos.get(lt)
    if not g: continue
    mid = (g["start"]+g["end"])/2
    y_gene = 0.5 if g["strand"]=="+" else -0.5
    y_label = LABEL_Y.get(lt, -1.7)
    x_label = mid + LABEL_NUDGE.get(lt, 0)
    ax1.annotate(label, xy=(mid, y_gene), xytext=(x_label, y_label),
                 ha="center", va="top", fontsize=7.5, color=col, fontweight="bold",
                 arrowprops=dict(arrowstyle="-", color=col, lw=0.6,
                                 connectionstyle="arc3,rad=0.1"), zorder=10)

ax1.set_xlim(-1500, pl25_L+1500)
ax1.set_ylim(-2.3, 2.3)
ax1.set_yticks([0.5, -0.5])
ax1.set_yticklabels(["+ strand", "− strand"], fontsize=9)
ax1.set_xticklabels([])
for s in ("top","right"): ax1.spines[s].set_visible(False)
ax1.set_title("A — Virgibacillus salarius PL25 contig_11 (74,025 bp) — 99 CDS, conjugative plasmid",
              fontsize=11, fontweight="bold", loc="left")

# ── Panel B: GC deviation ────────────────────────────────────────
ax2 = fig.add_subplot(gs[1], sharex=ax1)
xs = gc_pl25[:,0]; ys = gc_pl25[:,1]
deviation = ys - CHROM_GC

ax2.fill_between(xs, deviation, 0,
                 where=(deviation<=0), color="#8B1A1A", alpha=0.7, step=None)
ax2.fill_between(xs, deviation, 0,
                 where=(deviation>0), color="#1A6E3C", alpha=0.7, step=None)
ax2.axhline(0, color="#444", lw=0.7, ls="--")
mean_dev = deviation.mean()
ax2.axhline(mean_dev, color="black", lw=0.7, ls=":", alpha=0.6)
ax2.text(pl25_L*0.98, mean_dev, f"  mean Δ = {mean_dev:.1f}%", fontsize=7,
         va="center", ha="right", color="black")
ax2.set_xlim(-1500, pl25_L+1500)
ax2.set_ylim(-12, 4)
ax2.set_ylabel("GC – 36.81%\n(chrom mean)", fontsize=9)
ax2.set_xticklabels([])
for s in ("top","right"): ax2.spines[s].set_visible(False)
ax2.grid(axis="x", alpha=0.2)
ax2.set_title("B — GC deviation from chromosome (HGT signature)",
              fontsize=10, fontweight="bold", loc="left")

# ── Panel C: synteny ribbons ──────────────────────────────────
ax3 = fig.add_subplot(gs[2], sharex=ax1)

# Layout: PL25 top at y=1, O. contaminans bottom at y=0
ax3.set_ylim(-0.05, 1.05)
ax3.set_yticks([0, 1])
ax3.set_yticklabels(["O. contaminans\n(NZ_LDPV02000047.1)", "PL25 contig_11\n(this study)"], fontsize=8)

# Compute scaling factor: O. contaminans is 60,767 bp vs PL25 74,025 bp
# Show both on same x-axis using PL25 coords, scale O. contaminans linearly
scale = pl25_L / oc_L  # >1 because PL25 is bigger; we'll stretch oc to match

n_drawn = 0
for q, t, ident in alignments:
    if q not in pl25_pos or t not in oc_pos: continue
    qs, qe = pl25_pos[q]["start"], pl25_pos[q]["end"]
    ts, te = oc_pos[t]["start"]*scale, oc_pos[t]["end"]*scale
    # Colour by identity
    if ident >= 0.7: col = "#1F77B4"
    elif ident >= 0.5: col = "#5599DD"
    else: col = "#9BC9F0"
    alpha = 0.25 + 0.55*ident
    
    # Build polygon ribbon
    n = 30
    tlin = np.linspace(0, 1, n)
    x1 = qs*(1-tlin) + ts*tlin
    y1 = 1.0*(1-tlin) + 0.0*tlin
    x2 = qe*(1-tlin) + te*tlin
    y2 = y1
    pts_x = list(x1) + list(x2[::-1])
    pts_y = list(y1) + list(y2[::-1])
    ax3.fill(pts_x, pts_y, color=col, alpha=alpha, lw=0, zorder=2)
    n_drawn += 1

ax3.set_xticklabels([])
for s in ("top","right"): ax3.spines[s].set_visible(False)
ax3.set_title(f"C — Synteny: PL25 contig_11 ↔ O. contaminans (closest characterised relative)  •  {n_drawn} orthologous ribbons",
              fontsize=10, fontweight="bold", loc="left")

# ── Panel D: O. contaminans gene map ──────────────────────────
ax4 = fig.add_subplot(gs[3], sharex=ax1)
ax4.plot([0, pl25_L], [0, 0], color="#444", lw=1.0)
for g in oc_genes:
    y = 0.5 if g["strand"]=="+" else -0.5
    s_scaled = g["start"]*scale
    e_scaled = g["end"]*scale
    # Check if this protein is in the PL25 alignment set
    in_synteny = any(t == g["name"] for q, t, _ in alignments)
    col = "#444" if in_synteny else "#CCC"
    alpha = 0.85 if in_synteny else 0.4
    draw_gene(ax4, s_scaled, e_scaled, y, g["strand"], col, height=0.40,
              edgecolor="none", alpha=alpha)
ax4.set_xlim(-1500, pl25_L+1500)
ax4.set_ylim(-1.0, 1.0)
ax4.set_yticks([0.5,-0.5])
ax4.set_yticklabels(["+ strand","− strand"], fontsize=9)
ax4.set_xticklabels([])
for s in ("top","right"): ax4.spines[s].set_visible(False)
ax4.set_title("D — Ornithinibacillus contaminans contig (60,767 bp; rescaled to PL25 coordinates)",
              fontsize=10, fontweight="bold", loc="left")

# ── Panel E: position axis + summary box ───────────────────────────
ax5 = fig.add_subplot(gs[4], sharex=ax1)
ax5.set_xlim(-1500, pl25_L+1500)
ax5.set_ylim(0, 1)
ax5.set_yticks([])
ax5.set_xlabel("Position on contig_11 (bp)", fontsize=11)
for s in ("top","right","left"): ax5.spines[s].set_visible(False)

# Add a summary text box
summary = (
    "Virgibacillus salarius PL25 conjugative plasmid contig_11 (74,025 bp)\n"
    f"GC = 31.2% (chromosome 36.8%, Δ = −5.6%) — strong HGT signature\n"
    "geNomad plasmid score = 0.989  •  MOB-suite typing: novel (no canonical relaxase/MPF)\n"
    "M23 peptidase (PL25_00078) embedded in conjugation operon; closest characterised relative: O. contaminans (66.6% protein identity)"
)
ax5.text(0.5, 0.5, summary, transform=ax5.transAxes,
         ha="center", va="center", fontsize=9, color="#333",
         bbox=dict(facecolor="#FFF9E6", edgecolor="#C0A000", boxstyle="round,pad=0.4", lw=1))

# ── Title ─────────────────────────────────────────────────────────
fig.suptitle("Mobile genetic element architecture — V. salarius PL25 contig_11 conjugative plasmid",
             fontsize=13, fontweight="bold", y=0.995)

# ── Legend (mobility genes) ─────────────────────────────────────
legend_items = [
    mpatches.Patch(color="#9C27B0", label="Replication-Relaxation"),
    mpatches.Patch(color="#673AB7", label="RusA resolvase"),
    mpatches.Patch(color="#D7191C", label="M23 peptidase ★"),
    mpatches.Patch(color="#1976D2", label="VirB4-like ATPase"),
    mpatches.Patch(color="#42A5F5", label="ATPase / helicase"),
    mpatches.Patch(color="#5E3C99", label="T4SS-DNA transfer"),
    mpatches.Patch(color="#BBB",    label="Other plasmid CDS"),
    mpatches.Patch(color="#1F77B4", label="Synteny ≥ 70% id"),
    mpatches.Patch(color="#5599DD", label="Synteny 50–70%"),
    mpatches.Patch(color="#9BC9F0", label="Synteny 30–50%"),
    mpatches.Patch(color="#1A6E3C", label="GC > chrom"),
    mpatches.Patch(color="#8B1A1A", label="GC < chrom"),
]
fig.legend(handles=legend_items, loc="lower center", ncol=6,
           fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.01))

plt.tight_layout(rect=[0, 0.03, 1, 0.97])
fig.savefig(ROOT/"figures/07_mobile_element_integrated.pdf", bbox_inches="tight")
fig.savefig(ROOT/"figures/07_mobile_element_integrated.svg", bbox_inches="tight")
fig.savefig(ROOT/"figures/07_mobile_element_integrated.png", dpi=200, bbox_inches="tight")
print(f"✓ Saved 07_mobile_element_integrated.{{pdf,svg,png}}")
