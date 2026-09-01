#!/usr/bin/env python3
"""
Compact M23 tree with targeted keep-list (manuscript Fig 2 version).

Always kept uncollapsed:
  * PL25_M23
  * WP_066194568.1  ->  Cytobacillus horneckiae
  * WP_047980627.1  ->  Ornithinibacillus contaminans
  * REF_CwlP_clipped... (CwlP catalytic core)
  * The 5 REF M23 enzymes (Lysostaphin, LytM, ALE-1, EnpA, ZoocinA)

Phage M23 lysins are collapsed into a single tip 'Phage M23 lysins (n=14)'.

The remaining ~78 Bacillaceae homologues are collapsed into maximal
monophyletic wedges and labelled by dominant genus.

Inputs:
  results/07_M23_phylogeny/M23_iqtree_v4.contree
  results/07_M23_phylogeny/M23_seedset_v3.faa
Outputs:
  results/07_M23_phylogeny/M23_compact_targeted.nwk
  results/07_M23_phylogeny/M23_compact_targeted_summary.tsv
  figures/02e_M23_phylogeny_compact_targeted.{pdf,svg,png}
"""
import re
import copy
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from Bio import Phylo, SeqIO

ROOT = Path(__file__).resolve().parent.parent
PHY  = ROOT / "results" / "07_M23_phylogeny"
FIG  = ROOT / "figures"

CONTREE = PHY / "M23_iqtree_v4.contree"
SEEDSET = PHY / "M23_seedset_v3.faa"
OUT_NWK = PHY / "M23_compact_targeted.nwk"
OUT_TSV = PHY / "M23_compact_targeted_summary.tsv"

# ──────────────────────────────────────────────────────────────────
# Must-keep tips (these stay uncollapsed)
MUST_KEEP_PREFIX = ("PL25_M23", "REF_")
MUST_KEEP_EXACT = {
    "WP_066194568.1",   # Cytobacillus horneckiae
    "WP_047980627.1",   # Ornithinibacillus contaminans
}

# Category colours
COL = {
    "PL25":  "#D62728",
    "REF":   "#6A3D9A",
    "CwlP":  "#6A3D9A",
    "PHAGE": "#FF7F00",
    "BAC":   "#000000",       # bacillaceae wedges – plain black
    "BAC_HIGHLIGHT": "#1F77B4",  # the two preserved BAC tips
}
UFB_COL = [(95, "#1A7A1A", 30), (85, "#9A9A33", 22), (70, "#A0A0A0", 14)]

# ──────────────────────────────────────────────────────────────────
desc_map = {r.id: r.description for r in SeqIO.parse(str(SEEDSET), "fasta")}

def genus_of(name):
    if name.startswith("PHAGE_"):
        return "Phage"
    desc = desc_map.get(name, "")
    m = re.search(r"\[([A-Z][a-z]+)", desc)
    return m.group(1) if m else "Bacillaceae"

def is_must_keep(n):
    return n and (n.startswith(MUST_KEEP_PREFIX) or n in MUST_KEEP_EXACT)

def short_label(name):
    if name == "PL25_M23":
        return "PL25_M23 (this study)"
    if name == "REF_CwlP_clipped_SLT_M23_BACSU_SPbeta":
        return "CwlP (B. subtilis SP-β prophage)"
    if name.startswith("REF_"):
        return name.replace("REF_", "").split("_")[0]
    if name == "WP_066194568.1":
        return "WP_066194568.1 | Cytobacillus horneckiae"
    if name == "WP_047980627.1":
        return "WP_047980627.1 | Ornithinibacillus contaminans"
    return name

# ──────────────────────────────────────────────────────────────────
# Load + midpoint-root + ladderise
tree = Phylo.read(str(CONTREE), "newick")
tree.root_at_midpoint()
tree.ladderize()
print(f"loaded {len(tree.get_terminals())} tips; midpoint-rooted; ladderised")

# ──────────────────────────────────────────────────────────────────
# Step 1: collapse all 14 PHAGE_ tips into a single 'Phage M23 lysins (n=14)'
# placeholder attached at the MRCA of the original phage tips.
phage_tips = [t for t in tree.get_terminals() if t.name and t.name.startswith("PHAGE_")]
n_phage    = len(phage_tips)
if phage_tips:
    mrca = tree.common_ancestor(phage_tips)
    # Mean tip-to-mrca distance to give the synthetic terminal a sensible
    # branch length
    depths = tree.depths()
    mean_extra = float(np.mean([depths[t] for t in phage_tips]) - depths[mrca])

    # Prune each phage tip from the tree
    for t in list(phage_tips):
        try:
            tree.prune(t)
        except Exception:
            pass

    # The MRCA may have been pruned away if all its descendants were phages.
    # In that case the surviving parent of the MRCA position now becomes
    # the attachment point. Re-search for the MRCA's surviving location:
    surviving = set(tree.get_nonterminals())
    if mrca not in surviving:
        mrca = tree.root        # fallback — attach near the root

    # Create the synthetic terminal
    synthetic = Phylo.BaseTree.Clade(name=f"Phage M23 lysins (n={n_phage})",
                                     branch_length=max(0.05, mean_extra))
    mrca.clades.append(synthetic)
    print(f"  phages: pruned {n_phage} tips → 1 synthetic terminal at MRCA")

# Re-ladderise after surgery
tree.ladderize()

# ──────────────────────────────────────────────────────────────────
# Step 2: collapse remaining Bacillaceae homologues by monophyletic wedge
def descendant_leaf_names(c):
    return [t.name for t in c.get_terminals()]

# treat the synthetic phage terminal as must-keep so it survives wedge-collapse
def is_must_keep2(n):
    return is_must_keep(n) or (n and n.startswith("Phage M23 lysins"))

def can_collapse(c):
    return all(not is_must_keep2(n) for n in descendant_leaf_names(c))

wedges = []
def collapse_clade(c):
    leaves = descendant_leaf_names(c)
    gen = Counter(genus_of(n) for n in leaves)
    dom, dn = gen.most_common(1)[0]
    n = len(leaves)
    purity = dn / n
    if purity >= 0.75:
        lab = f"{dom} clade (n={n})"
    elif purity >= 0.50:
        top2 = gen.most_common(2)
        lab = f"{top2[0][0]} / {top2[1][0]} clade (n={n})"
    else:
        lab = f"Bacillaceae mixed clade (n={n})"
    # Use mean tip depth as the wedge's branch length so the layout still scales
    depths = tree.depths()
    extra = float(np.mean([depths[t] for t in c.get_terminals()]) - depths[c])
    c.branch_length = (c.branch_length or 0.0) + max(extra, 0.0)
    c.clades = []
    c.name = lab
    c.confidence = None
    wedges.append({"label": lab, "n_leaves": n, "dominant_genus": dom,
                   "dominant_fraction": purity, "leaves": ";".join(leaves)})

def walk(c):
    if c.is_terminal():
        return
    if can_collapse(c) and len(descendant_leaf_names(c)) >= 2:
        collapse_clade(c)
        return
    for ch in list(c.clades):
        walk(ch)

walk(tree.root)
tree.ladderize()

# Rewrite tip names: short labels for must-keep; the rest already have wedge labels
for t in tree.get_terminals():
    if t.name and is_must_keep(t.name):
        t.name = short_label(t.name)

# Save Newick + wedge summary
Phylo.write(tree, str(OUT_NWK), "newick")
with OUT_TSV.open("w") as fh:
    fh.write("label\tn_leaves\tdominant_genus\tdominant_fraction\tleaves\n")
    for w in wedges:
        fh.write(f"{w['label']}\t{w['n_leaves']}\t{w['dominant_genus']}\t"
                 f"{w['dominant_fraction']:.2f}\t{w['leaves']}\n")
n_remaining = len(tree.get_terminals())
print(f"  Bacillaceae wedges created: {len(wedges)}")
print(f"  total tips remaining      : {n_remaining}")
print(f"  newick                    : {OUT_NWK}")
print(f"  wedge summary             : {OUT_TSV}")

# ──────────────────────────────────────────────────────────────────
# Render — same visual language as scripts/25
def category(name):
    if not name: return "BAC"
    if name.startswith("PL25_M23"):        return "PL25"
    if name.startswith("CwlP"):            return "CwlP"
    if name in {"Lysostaphin","LytM","ALE-1","EnpA","ZoocinA"}: return "REF"
    if name.startswith("Phage M23 lysins"):return "PHAGE"
    if name.startswith("WP_"):             return "BAC_HIGHLIGHT"
    return "BAC"

def tip_color(name):
    cat = category(name)
    return COL.get(cat, "#000000")

def tip_marker(name):
    cat = category(name)
    if cat == "PL25":          return ("*", 200, COL["PL25"], "black", 1.0)
    if cat == "CwlP":          return ("s", 60,  COL["CwlP"], "black", 0.5)
    if cat == "REF":           return ("s", 55,  COL["REF"],  "black", 0.5)
    if cat == "PHAGE":         return ("D", 70,  COL["PHAGE"],"black", 0.5)
    if cat == "BAC_HIGHLIGHT": return ("o", 30,  COL["BAC_HIGHLIGHT"], "black", 0.4)
    return (None, 0, None, None, 0)

# Coordinates
leaves = tree.get_terminals()
N = len(leaves)
depths = tree.depths()
y_of = {tip: N - i for i, tip in enumerate(leaves)}
def assign_internal(c):
    for ch in c.clades:
        if ch not in y_of:
            assign_internal(ch)
    y_of[c] = sum(y_of[ch] for ch in c.clades) / len(c.clades)
assign_internal(tree.root)

MAX_X = max(depths.values())
fig_h = max(7.0, N * 0.32)
fig, ax = plt.subplots(figsize=(8.0, fig_h))

# Draw branches
for c in tree.find_clades():
    parent = None
    for p in tree.find_clades():
        if p is c: continue
        if c in p.clades:
            parent = p; break
    if parent is None: continue
    x0, x1 = depths[parent], depths[c]
    y0, y1 = y_of[parent], y_of[c]
    ax.plot([x0, x0], [y0, y1], color="#444", lw=0.7, zorder=2)
    ax.plot([x0, x1], [y1, y1], color="#444", lw=0.7, zorder=2)

# Draw tips
TIP_PAD = MAX_X * 0.02

def italicize_clade_label(label):
    """Italicise genus/family names in collapsed-clade labels via mathtext.
    Handles:
      'X / Y clade (n=N)'      → italic X, italic Y
      'X clade (n=N)'          → italic X
      'Bacillaceae mixed …'    → italic Bacillaceae
    Non-matching strings returned unchanged.
    """
    m = re.match(r"^([A-Z][a-z]+) / ([A-Z][a-z]+) clade \(n=(\d+)\)$", label)
    if m:
        g1, g2, n = m.groups()
        return rf"$\it{{{g1}}}$ / $\it{{{g2}}}$ clade (n={n})"
    m = re.match(r"^([A-Z][a-z]+) clade \(n=(\d+)\)$", label)
    if m:
        g, n = m.groups()
        return rf"$\it{{{g}}}$ clade (n={n})"
    m = re.match(r"^(Bacillaceae) mixed clade \(n=(\d+)\)$", label)
    if m:
        fam, n = m.groups()
        return rf"$\it{{{fam}}}$ mixed clade (n={n})"
    return label

for tip in leaves:
    x = depths[tip]
    y = y_of[tip]
    cat = category(tip.name)
    marker, size, fc, ec, lw = tip_marker(tip.name)
    if marker is not None:
        ax.scatter([x + TIP_PAD*0.5], [y], s=size, marker=marker,
                   facecolor=fc, edgecolor=ec, linewidths=lw, zorder=5)

    label = tip.name
    if cat == "PL25":
        ax.text(x + TIP_PAD*2.0, y, label, fontsize=10, color=COL["PL25"],
                va="center", fontweight="bold", zorder=6)
    elif cat in ("CwlP", "REF"):
        ax.text(x + TIP_PAD*2.0, y, label, fontsize=9, color=COL["REF"],
                va="center", fontweight="bold", zorder=6)
    elif cat == "PHAGE":
        ax.text(x + TIP_PAD*2.0, y, label, fontsize=9, color=COL["PHAGE"],
                va="center", zorder=6)
    elif cat == "BAC_HIGHLIGHT":
        # WP_066194568.1 | Cytobacillus horneckiae   etc.
        ax.text(x + TIP_PAD*1.6, y, label, fontsize=9, color="#000000",
                va="center", zorder=6)
    else:
        ax.text(x + TIP_PAD*0.5, y, italicize_clade_label(label),
                fontsize=9, color="#000000", va="center", zorder=4)

# UFBoot circles
for c in tree.find_clades(terminal=False):
    ufb = None
    if c.confidence is not None:
        try:    ufb = float(c.confidence)
        except: ufb = None
    if ufb is None and c.name is not None:
        try:    ufb = float(c.name)
        except: ufb = None
    if ufb is None or ufb < 70: continue
    for thresh, col, sz in UFB_COL:
        if ufb >= thresh:
            ax.scatter([depths[c]], [y_of[c]],
                       s=sz, marker="o", facecolor=col,
                       edgecolor=col, linewidths=0.4, zorder=6)
            break

# Scale bar
scale_len = round(MAX_X / 5, 1)
if scale_len <= 0: scale_len = MAX_X / 5
sx = MAX_X * 0.02
sy = -1.0
ax.plot([sx, sx + scale_len], [sy, sy], "k-", lw=1.2)
ax.text(sx + scale_len/2, sy - 0.5, f"{scale_len} substitutions/site",
        ha="center", va="top", fontsize=8)

ax.set_xlim(-MAX_X*0.02, MAX_X * 1.55)
ax.set_ylim(-1.5, N + 1)
ax.axis("off")

# Legend
handles = [
    Line2D([0],[0], marker="*", linestyle="", color=COL["PL25"],
           markersize=14, markeredgecolor="black",
           label="PL25_M23 (this study)"),
    Line2D([0],[0], marker="s", linestyle="", color=COL["REF"],
           markersize=8, markeredgecolor="black",
           label="Reference M23 enzymes / CwlP"),
    Line2D([0],[0], marker="D", linestyle="", color=COL["PHAGE"],
           markersize=8, markeredgecolor="black",
           label="Phage M23 lysins (collapsed)"),
    Line2D([0],[0], marker="o", linestyle="", color=COL["BAC_HIGHLIGHT"],
           markersize=7, markeredgecolor="black",
           label="Key Bacillaceae neighbours"),
    Line2D([0],[0], marker="o", linestyle="", color="#1A7A1A",
           markersize=9, label="UFBoot ≥ 95"),
    Line2D([0],[0], marker="o", linestyle="", color="#9A9A33",
           markersize=7, label="UFBoot 85–94"),
    Line2D([0],[0], marker="o", linestyle="", color="#A0A0A0",
           markersize=6, label="UFBoot 70–84"),
]
ax.legend(handles=handles, loc="lower right", fontsize=8.5,
          frameon=True, framealpha=0.95, handletextpad=0.6, borderpad=0.7)

fig.tight_layout()
stem = FIG / "02e_M23_phylogeny_compact_targeted"
for ext in ("pdf","svg","png"):
    fig.savefig(f"{stem}.{ext}", bbox_inches="tight", dpi=300)
print(f"  figure                    : {stem}.{{pdf,svg,png}}")
