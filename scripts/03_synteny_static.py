#!/usr/bin/env python3
"""Static synteny figure: PL25 contig_11 vs 3 close NCBI relatives.
Reproduces the style of 12_static_synteny.py from the prophage paper.
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO

ROOT = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
SYNT = ROOT/"results/14_synteny"

tracks = [
    ("PL25_contig11",          SYNT/"PL25_contig11.gbk",          "$\\it{V.\\ salarius}$ PL25 (this study)",                 "#D7191C"),
    ("NZ_LDPV02000047.1",      SYNT/"NZ_LDPV02000047.1.gbk",      "$\\it{Ornithinibacillus\\ contaminans}$ (closest)",       "#1F77B4"),
    ("NZ_JARSFA010000044.1",   SYNT/"NZ_JARSFA010000044.1.gbk",   "$\\it{Cytobacillus\\ horneckiae}$",                       "#2CA02C"),
    ("NZ_JBCITH010000002.1",   SYNT/"NZ_JBCITH010000002.1.gbk",   "$\\it{Caldifermentibacillus\\ hisashii}$",                "#FF7F0E"),
]

# Parse GBK files
def parse_gbk(path):
    rec = next(SeqIO.parse(str(path),"genbank"))
    genes = []
    for f in rec.features:
        if f.type != "CDS": continue
        q = f.qualifiers
        lt = q.get("protein_id", q.get("locus_tag",["?"]))[0]
        if not lt or lt == "?":
            lt = q.get("protein_id",["?"])[0]
        prod = q.get("product",["hypothetical"])[0]
        try:
            s = int(f.location.start)
            e = int(f.location.end)
        except: continue
        strand = "+" if f.location.strand == 1 else "-"
        genes.append({"name": lt, "start": s, "end": e, "strand": strand, "product": prod})
    return genes, len(rec.seq)

track_genes = {}
track_len = {}
for tid, path, _, _ in tracks:
    g, L = parse_gbk(path)
    track_genes[tid] = g
    track_len[tid] = L
    print(f"  {tid}: {L:,} bp, {len(g)} CDS")

max_len = max(track_len.values())

# Categorise key genes for colouring
KEY_GENES = {
    "PL25_00076": ("Replication-Relaxation", "#9C27B0"),
    "PL25_00077": ("RusA resolvase",         "#673AB7"),
    "PL25_00078": ("M23 peptidase",          "#D7191C"),
    "PL25_00082": ("VirB4-like ATPase",      "#1F77B4"),
    "PL25_00094": ("ATPase",                 "#42A5F5"),
    "PL25_00148": ("T4SS-DNA transfer",      "#5E3C99"),
}

# Read clinker alignments  
alignments = []
in_section = False
cur_q, cur_t = None, None
with open(SYNT/"PL25_vs_relatives_alignments.csv") as f:
    for raw in f:
        line = raw.rstrip()
        if not line.strip():
            in_section = False; continue
        if line.startswith("---"):
            continue
        if " vs " in line and not line.startswith("Query"):
            parts = line.split(" vs ")
            cur_q = parts[0].strip(); cur_t = parts[1].strip()
            in_section = False; continue
        parts = line.split()
        if len(parts) < 4: continue
        if parts[0] == "Query" and parts[1] == "Target":
            in_section = True; continue
        if in_section:
            try:
                q, t, ident = parts[0], parts[1], float(parts[2])
                alignments.append((cur_q, q, cur_t, t, ident))
            except: pass

print(f"  Alignments: {len(alignments)}")

# Build gene-position lookup
pos_map = {}
for tid, _, _, _ in tracks:
    for g in track_genes[tid]:
        pos_map[g["name"]] = (tid, g["start"], g["end"])

# Draw
def draw_gene(ax, s, e, y, strand, col, height=0.45):
    arrow_w = min(800, (e-s)*0.35)
    body_w = (e-s) - arrow_w
    if strand == "+":
        x = [s, s+body_w, s+body_w, e, s+body_w, s+body_w, s, s]
        ydir = [y-height/2, y-height/2, y-height, y, y+height, y+height/2, y+height/2, y-height/2]
    else:
        x = [e, s+arrow_w, s+arrow_w, s, s+arrow_w, s+arrow_w, e, e]
        ydir = [y-height/2, y-height/2, y-height, y, y+height, y+height/2, y+height/2, y-height/2]
    ax.fill(x, ydir, color=col, lw=0.4, edgecolor="black", zorder=4)

def draw_ribbon(ax, qs, qe, yq, ts, te, yt, col, alpha):
    # Curved ribbon between (qs,qe,yq) and (ts,te,yt)
    import numpy as np
    n = 50
    t = np.linspace(0, 1, n)
    mid_y = (yq + yt) / 2
    x1 = qs*(1-t) + ts*t
    y1 = yq*(1-t) + yt*t
    x2 = qe*(1-t) + te*t
    y2 = yq*(1-t) + yt*t
    ax.fill_betweenx(np.concatenate([y1, y2[::-1]]),
                     np.concatenate([x1, x2[::-1]]),
                     np.concatenate([x1, x2[::-1]]),
                     facecolor=col, alpha=alpha, lw=0, zorder=2)
    # Simpler: just polygon  
    pts_x = list(x1) + list(x2[::-1])
    pts_y = list(y1) + list(y2[::-1])
    ax.fill(pts_x, pts_y, color=col, alpha=alpha, lw=0, zorder=2)

n = len(tracks)
fig, ax = plt.subplots(figsize=(16, 1.6*n + 1.5))
track_y = {}
for i, (tid, _, label, col) in enumerate(tracks):
    y = (n-1-i) * 1.6
    track_y[tid] = y
    ax.plot([0, track_len[tid]], [y, y], color="#444", lw=1.0, zorder=1)
    ax.text(-max_len*0.012, y, label,
            ha="right", va="center", fontsize=8.5,
            fontweight="bold" if i==0 else "normal",
            color=col if i==0 else "#222")
    for g in track_genes[tid]:
        if g["name"] in KEY_GENES:
            gc = KEY_GENES[g["name"]][1]
        else:
            gc = "#999"
        draw_gene(ax, g["start"], g["end"], y, g["strand"], gc, height=0.40)

# --- Build orthology graph so we can trace any gene back to a PL25 key gene -----
# Two-pass: first collect undirected pairs, then BFS from each KEY_GENES seed.
from collections import defaultdict, deque
adj = defaultdict(set)
for cur_q, q, cur_t, t, ident in alignments:
    if q in pos_map and t in pos_map:
        adj[q].add(t); adj[t].add(q)
linked_to_key = set(KEY_GENES.keys())
queue = deque(KEY_GENES.keys())
while queue:
    cur = queue.popleft()
    for nb in adj[cur]:
        if nb not in linked_to_key:
            linked_to_key.add(nb); queue.append(nb)

# Ribbons — bright for any ribbon whose endpoints are in the key-linked set,
#          dimmed grey for everything else (kept for context but de-emphasised).
n_rib_key = 0; n_rib_bg = 0
for cur_q, q, cur_t, t, ident in alignments:
    if q not in pos_map or t not in pos_map: continue
    tq, qs, qe = pos_map[q]
    tt, ts, te = pos_map[t]
    if tq not in track_y or tt not in track_y: continue
    if abs([x[0] for x in tracks].index(tq) - [x[0] for x in tracks].index(tt)) != 1: continue
    yq = track_y[tq]; yt = track_y[tt]
    is_key = (q in linked_to_key) and (t in linked_to_key)
    if is_key:
        # propagate colour from the closest PL25 key gene if we can find it on this track
        col = KEY_GENES.get(q, KEY_GENES.get(t, (None, "#D7191C")))[1]
        alpha = 0.35 + 0.45*float(ident)
        n_rib_key += 1
    else:
        col = "#BBBBBB"
        alpha = 0.08
        n_rib_bg += 1
    draw_ribbon(ax, qs, qe, yq, ts, te, yt, col, alpha)
print(f"  Ribbons drawn: key={n_rib_key}, background={n_rib_bg}")

# (Per-gene labels on the PL25 track removed; gene identities are now shown in the legend.)

ax.set_xlim(-max_len*0.20, max_len*1.02)
ax.set_ylim(-1, n*1.6 + 0.5)
ax.set_yticks([])
ax.set_xlabel("Genomic position (bp)")
for s in ("top","right","left"): ax.spines[s].set_visible(False)
ax.grid(axis="x", alpha=0.25)

ax.set_title("Synteny — V. salarius PL25 contig_11 vs closest halophile-Bacillaceae relatives\n"
             "clinker (identity ≥ 0.30) • 4 sequences • PL25 conjugation operon shows strong synteny with Ornithinibacillus contaminans",
             fontsize=12, fontweight="bold", pad=12)

# Legend — locus tag + functional descriptor side-by-side
legend_items = [mpatches.Patch(color=col, label=f"{lt}  {desc}")
                for lt, (desc, col) in KEY_GENES.items()]
legend_items.append(mpatches.Patch(color="#999", label="Other CDS"))
ax.legend(handles=legend_items, loc="upper right", bbox_to_anchor=(1.0, -0.06),
          ncol=4, fontsize=8.5, frameon=True, framealpha=0.9, title="Key PL25 genes")

# Footnote
fig.text(0.5, 0.01,
    "Independent Prokka annotation of PL25 contig_11 vs three close halotolerant Bacillaceae contigs from NCBI nr. "
    "Each NCBI accession (full GBK record) was selected as the contig containing the M23 protein with highest sequence similarity to PL25_M23. "
    "Strong block synteny over the conjugation operon (left side) confirms shared evolutionary origin.",
    ha="center", fontsize=7.5, style="italic", color="#444", wrap=True)

plt.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(ROOT/"figures/06_synteny_contig11_vs_relatives.pdf", bbox_inches="tight")
fig.savefig(ROOT/"figures/06_synteny_contig11_vs_relatives.svg", bbox_inches="tight")
fig.savefig(ROOT/"figures/06_synteny_contig11_vs_relatives.png", dpi=200, bbox_inches="tight")
print(f"✓ Saved 06_synteny_contig11_vs_relatives.{{pdf,svg,png}}")
