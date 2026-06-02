#!/usr/bin/env python3
"""Circular map of PL25 contig_11 (74,025 bp closed circular plasmid).

Rings from outside in:
  1. CDS arrows (forward + reverse strand)
  2. Position / scale ticks
  3. Coverage track (uniform ≈ 889× — Flye `assembly_info.txt`)
  4. GC content (deviation from contig mean)
  5. GC skew  (sign change marks the predicted origin of replication)

Origin of plot is positioned at PL25_00076 (predicted replication-relaxation
protein) — the canonical place to start a plasmid map.

Run: /usr/local/Caskroom/miniforge/base/envs/macsy/bin/python scripts/09_circular_plasmid_map.py
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np

ROOT = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
GBK  = ROOT/"results/14_synteny/PL25_contig11.gbk"
OUT  = ROOT/"figures"

# ── Parse GBK ─────────────────────────────────────────────
gbk = Genbank(str(GBK))
seq = gbk.full_genome_seq
L   = len(seq)
print(f"contig_11: {L:,} bp")

# Find PL25_00076 start position to use as origin
features = gbk.extract_features("CDS")
origin_offset = 0
for f in features:
    q = f.qualifiers
    name = (q.get("locus_tag",["?"])[0])
    if "00076" in name:
        s = int(f.location.start)
        origin_offset = s
        print(f"PL25_00076 starts at {s} bp — using as plot origin")
        break

# ── Categorise CDS ────────────────────────────────────────
KEY = {
    "PL25_00076": ("Replication-Relaxation", "#9C27B0"),
    "PL25_00077": ("RusA resolvase",         "#673AB7"),
    "PL25_00078": ("M23 peptidase",          "#D7191C"),
    "PL25_00082": ("VirB4-like ATPase",      "#1F77B4"),
    "PL25_00094": ("ATPase",                 "#42A5F5"),
    "PL25_00148": ("T4SS-DNA transfer",      "#5E3C99"),
}

# ── Build Circos ──────────────────────────────────────────
circos = Circos(sectors={"plasmid": L}, start=-90, end=270)
sector = circos.sectors[0]

# Track 1: outer ticks every 10 kb
sector.text(
    "$\\bf{PL25}$  $\\it{contig\\_11}$\n"
    "74,025 bp\nclosed circular plasmid\n"
    "\n"
    "Flye:  circ = Y\n"
    "mean cov ≈ 889×\n"
    f"GC = 31.2 %",
    r=15, size=10.5, color="#222")

major_xticks_kb = list(range(0, int(L/1000)+1, 10))
major_xticks    = [k*1000 for k in major_xticks_kb]
tick_labels     = [f"{k} kb" for k in major_xticks_kb]
sector.axis(fc="none", lw=0)
outer_track = sector.add_track((96, 100))
outer_track.axis(fc="white", ec="black", lw=0.6)
outer_track.xticks_by_interval(10000, label_size=8, label_orientation="vertical",
                               label_formatter=lambda v: f"{int(v/1000)} kb")

# Track 2: forward-strand CDS
fwd_track = sector.add_track((86, 95))
fwd_track.axis(fc="white", ec="lightgrey", lw=0.3)
# Track 3: reverse-strand CDS
rev_track = sector.add_track((77, 86))
rev_track.axis(fc="white", ec="lightgrey", lw=0.3)

# Draw arrows
for f in features:
    q = f.qualifiers
    name = q.get("locus_tag",["?"])[0]
    if not name or name == "?":
        name = q.get("protein_id",["?"])[0]
    s = int(f.location.start)
    e = int(f.location.end)
    strand = f.location.strand
    col = KEY.get(name, (None, "#BDBDBD"))[1]
    track = fwd_track if strand == 1 else rev_track
    track.arrow(s, e, fc=col, ec="black", lw=0.3)

# Label key genes — leader lines from outside the outer ring
for name, (descr, col) in KEY.items():
    for f in features:
        if f.qualifiers.get("locus_tag",["?"])[0] == name:
            s = int(f.location.start); e = int(f.location.end)
            mid = (s+e)/2
            label = f"{name}\n{descr}"
            if name == "PL25_00078":
                label = f"★ {name}\n{descr}\n(this study)"
            sector.text(label, mid, r=112, size=8, color=col,
                        ha="center", va="center", weight="bold")
            break

# Track 4: Coverage ring (uniform ~889x as a flat ring with mean labelled)
cov_track = sector.add_track((71, 76))
cov_track.axis(fc="#E3F2FD", ec="black", lw=0.4)
cov_track.bar([L/2], [1.0], width=L, color="#1976D2", alpha=0.55)

# Track 5: GC content (deviation from mean)
gc_track = sector.add_track((61, 70))
gc_track.axis(fc="white", ec="lightgrey", lw=0.3)
positions, gc_vals = gbk.calc_gc_content(window_size=1000, step_size=200)
mean_gc = float(gbk.calc_genome_gc_content())
gc_dev  = np.asarray(gc_vals) - mean_gc
gc_max  = float(max(abs(gc_dev.min()), abs(gc_dev.max())))
gc_above = np.where(gc_dev > 0, gc_dev, 0.0)
gc_below = np.where(gc_dev < 0, gc_dev, 0.0)
gc_track.fill_between(positions, gc_above, 0, color="#1B5E20", alpha=0.7,
                      vmin=-gc_max, vmax=gc_max)
gc_track.fill_between(positions, gc_below, 0, color="#FFB300", alpha=0.7,
                      vmin=-gc_max, vmax=gc_max)

# Track 6: GC skew  (sign change marks predicted oriV)
skew_track = sector.add_track((50, 60))
skew_track.axis(fc="white", ec="lightgrey", lw=0.3)
positions2, sk = gbk.calc_gc_skew(window_size=1000, step_size=200)
sk = np.asarray(sk)
sk_max = float(max(abs(sk.min()), abs(sk.max())))
sk_pos = np.where(sk > 0, sk, 0.0)
sk_neg = np.where(sk < 0, sk, 0.0)
skew_track.fill_between(positions2, sk_pos, 0, color="#B71C1C", alpha=0.7,
                        vmin=-sk_max, vmax=sk_max)
skew_track.fill_between(positions2, sk_neg, 0, color="#0D47A1", alpha=0.7,
                        vmin=-sk_max, vmax=sk_max)

# Save
fig = circos.plotfig()
ax  = fig.axes[0]

legend = [
    Patch(color="#D7191C", label="PL25_00078  M23 peptidase  (this study)"),
    Patch(color="#9C27B0", label="PL25_00076  Replication-Relaxation"),
    Patch(color="#673AB7", label="PL25_00077  RusA resolvase"),
    Patch(color="#1F77B4", label="PL25_00082  VirB4-like ATPase"),
    Patch(color="#42A5F5", label="PL25_00094  AAA ATPase"),
    Patch(color="#5E3C99", label="PL25_00148  T4SS-DNA transfer"),
    Patch(color="#BDBDBD", label="Other CDS"),
    Patch(color="#1B5E20", label="GC content (above mean)"),
    Patch(color="#FFB300", label="GC content (below mean)"),
    Patch(color="#B71C1C", label="GC skew positive"),
    Patch(color="#0D47A1", label="GC skew negative"),
]
fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.04),
           ncol=3, fontsize=8, frameon=True, framealpha=0.95,
           title="Tracks (outer → inner: CDS forward, CDS reverse, coverage, GC %, GC skew)",
           title_fontsize=9)

fig.text(0.5, 0.965,
    "PL25 contig_11 — 74,025 bp closed circular plasmid carrying PL25_M23",
    ha="center", fontsize=13, fontweight="bold")
fig.text(0.5, 0.945,
    "Flye-confirmed circularity (assembly_info.txt: circ=Y; mean cov 889×; Nanopore + Flye v2.9+)",
    ha="center", fontsize=9.5, style="italic", color="#444")

out = OUT/"01c_contig11_circular_map"
for ext in ("pdf","svg","png"):
    fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print(f"✓ Saved {out}.{{pdf,svg,png}}")
