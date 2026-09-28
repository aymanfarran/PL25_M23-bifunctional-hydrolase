#!/usr/bin/env python3
"""Render main-text Figure 3 (2-panel): focused phylogeny + replicon-context summary.

Inputs:
  results/22_closed_genome_search/closed_M23_iqtree.treefile   (produced by IQ-TREE)
Outputs:
  figures/F3_M23_plasmid_lineage.{pdf,svg,png}

Pipeline upstream (produces the IQ-TREE input):
  python scripts/22_closed_genome_search.py
  mafft --auto --thread 1 closed_M23_with_outgroup.faa > closed_M23_aln.fasta
  trimal -in closed_M23_aln.fasta -out closed_M23_aln_trim.fasta -gappyout
  iqtree3 -s closed_M23_aln_trim.fasta -m MFP -B 1000 -alrt 1000 \\
          -T 4 -pre closed_M23_iqtree

MAFFT is run single-threaded on purpose: its multithreaded iterative refinement
visits the guide tree in a nondeterministic order, and repeated --auto runs of
these 16 sequences gave alignments of 1069, 1071 and 1072 columns. With
--thread 1 the alignment, and so the tree, is identical on every run.

Run:
  /usr/local/Caskroom/miniforge/base/envs/macsy/bin/python \\
      scripts/22_render_F3_plasmid_lineage.py
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from Bio import Phylo
from Bio.Phylo.BaseTree import BranchColor
from collections import Counter

ROOT    = Path(__file__).resolve().parent.parent
TREE    = ROOT / "results/22_closed_genome_search/closed_M23_iqtree.treefile"
OUTBASE = ROOT / "figures/F3_M23_plasmid_lineage"

# ──────────────────────────────────────────────────────────
# Tree
# ──────────────────────────────────────────────────────────
tree = Phylo.read(str(TREE), "newick")
for leaf in tree.get_terminals():
    if leaf.name.startswith("Lysostaphin"):
        tree.root_with_outgroup(leaf); break

def name_to_color(name: str) -> str:
    if "PL25_M23"    in name: return "#D7191C"   # red — this study
    if "Lysostaphin" in name: return "#000000"   # black — outgroup
    if "|plasmid_"   in name: return "#1F77B4"   # blue — plasmid-encoded
    if "|chromosome" in name: return "#7F7F7F"   # grey — chromosomal
    return "#000000"

def hex_to_bc(h: str) -> BranchColor:
    h = h.lstrip("#"); return BranchColor(int(h[:2],16), int(h[2:4],16), int(h[4:],16))

for leaf in tree.get_terminals():
    leaf.color = hex_to_bc(name_to_color(leaf.name))

def propagate(c):
    if c.is_terminal(): return c.color
    kids = [propagate(k) for k in c.clades]
    if all((k.red,k.green,k.blue) == (kids[0].red,kids[0].green,kids[0].blue) for k in kids):
        c.color = kids[0]
    else:
        c.color = hex_to_bc("#222222")     # mixed
    return c.color
propagate(tree.root)

def parse_ufboot(clade):
    s = None
    if clade.confidence is not None: s = str(clade.confidence)
    elif clade.name: s = clade.name
    if not s or "/" not in s: return None
    try: return float(s.split("/")[1])
    except: return None

def tip_label(name: str) -> str:
    parts = name.split("|")
    raw = parts[0].replace("_", " "); words = raw.split()
    org = (words[0][0] + ". " + " ".join(words[1:])) if len(words) >= 2 else raw
    ctx = parts[1].replace("plasmid_", "").replace("chromosome", "chr")
    if "PL25_M23" in name:    return "★ PL25_M23  (pPL25-M23, this study)"
    if "Lysostaphin" in name: return "Lysostaphin  [outgroup]"
    return f"{org}   [{ctx}]"

pretty_to_orig = {tip_label(l.name): l.name for l in tree.get_terminals()}
def label_color(lbl):
    orig = pretty_to_orig.get(lbl)
    return name_to_color(orig) if orig else "#000"

# Y / X positions (same algorithm as Phylo.draw)
y_pos, x_pos = {}, {}
def assign_y(clade, counter):
    if clade.is_terminal():
        counter[0] += 1; y_pos[id(clade)] = counter[0]; return counter[0]
    ys = [assign_y(c, counter) for c in clade.clades]
    y_pos[id(clade)] = (min(ys)+max(ys))/2.0; return y_pos[id(clade)]
assign_y(tree.root, [0])
N = len(tree.get_terminals())
y_pos_flip = {k: (N - v + 1) for k, v in y_pos.items()}
def assign_x(clade, x):
    x_pos[id(clade)] = x
    for c in clade.clades: assign_x(c, x + (c.branch_length or 0))
assign_x(tree.root, 0)

# ──────────────────────────────────────────────────────────
# Panel-b summary
# ──────────────────────────────────────────────────────────
# Read from the census rather than restating it here, so the panel cannot drift
# away from the table and the tree if the search is ever rerun.
CENSUS = ROOT / "results/22_closed_genome_search/closed_M23_homologues.tsv"
TRIMMED = ROOT / "results/22_closed_genome_search/closed_M23_aln_trim.fasta"
with CENSUS.open() as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
hits = [(r["host"].split()[0], r["replicon"]) for r in rows]

# Alignment width quoted in the caption, measured rather than remembered.
N_COLS = len("".join(
    L.strip() for L in TRIMMED.read_text().split(">")[1].splitlines()[1:]))

plasmid = Counter(); chrom = Counter()
for g,r in hits: (plasmid if r=="plasmid" else chrom)[g] += 1
genera = sorted(set(plasmid)|set(chrom), key=lambda g: -(plasmid[g]+chrom[g]))

# ──────────────────────────────────────────────────────────
# Figure
# ──────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 7.2))
gs = GridSpec(1, 2, width_ratios=[2.4, 1], wspace=0.32)
ax_tree = fig.add_subplot(gs[0]); ax_bar = fig.add_subplot(gs[1])

Phylo.draw(tree, axes=ax_tree, do_show=False,
           label_func=lambda c: tip_label(c.name) if c.is_terminal() else "",
           label_colors=label_color, show_confidence=False)
for clade in tree.get_nonterminals():
    ub = parse_ufboot(clade)
    if ub is None: continue
    ax_tree.text(x_pos[id(clade)], y_pos_flip[id(clade)],
                 f"{int(ub)}", fontsize=7.5, ha="right", va="center",
                 color="#222",
                 bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
# Tip labels carry genus, species, strain and replicon, so the drawing area
# needs room to the right of the deepest tip or they run into panel b.
ax_tree.set_xlim(right=max(x_pos.values()) * 1.95)
ax_tree.set_title(f"a   Focused phylogeny of PL25_M23 and {len(hits)} closed-genome\n     halotolerant Bacillaceae homologues",
                  loc="left", fontsize=12, fontweight="bold", pad=6)
ax_tree.set_xlabel("Substitutions per site", fontsize=10); ax_tree.set_ylabel("")
ax_tree.legend(handles=[
    Patch(color="#D7191C", label="★ PL25_M23 (this study)"),
    Patch(color="#1F77B4", label="Plasmid-encoded (n = 11)"),
    Patch(color="#7F7F7F", label="Chromosomal (n = 3)"),
    Patch(color="#000000", label="Lysostaphin (outgroup)"),
], loc="lower right", fontsize=8, frameon=True, framealpha=0.95)

ys = range(len(genera))
plas_vals = [plasmid[g] for g in genera]; chr_vals = [chrom[g] for g in genera]
ax_bar.barh(ys, plas_vals, color="#1F77B4", edgecolor="#0F4F8A", label="plasmid")
ax_bar.barh(ys, chr_vals,  color="#7F7F7F", edgecolor="#404040",
            left=plas_vals, label="chromosome")
ax_bar.set_yticks(list(ys))
ax_bar.set_yticklabels([f"$\\it{{{g}}}$" for g in genera], fontsize=10)
ax_bar.invert_yaxis()
ax_bar.set_xlabel("Number of closed-genome homologues", fontsize=10)
ax_bar.set_xlim(0, max(p+c for p,c in zip(plas_vals, chr_vals)) + 0.5)
for i,(p,c) in enumerate(zip(plas_vals, chr_vals)):
    if p+c: ax_bar.text(p+c+0.1, i, f" {p+c}", va="center", fontsize=9, color="#222")
n_plasmid = sum(plasmid.values()); n_hits = len(hits)
ax_bar.set_title(f"b   Replicon context across genera\n"
                 f"     ({n_plasmid} of {n_hits}, {100*n_plasmid/n_hits:.0f}% plasmid-encoded)",
                 loc="left", fontsize=12, fontweight="bold", pad=6)
ax_bar.legend(loc="lower right", fontsize=9, frameon=True, framealpha=0.95)
ax_bar.spines["top"].set_visible(False); ax_bar.spines["right"].set_visible(False)
ax_bar.grid(axis="x", alpha=0.3)

fig.text(0.5, 0.01,
    f"(a) MAFFT v7.526 (--auto) → trimAl v1.5.rev1 (-gappyout, {N_COLS} cols) → "
    "IQ-TREE v3.1.1 (WAG+G4, 1000 UFBoot / 1000 SH-aLRT). Numbers at internal nodes are UFBoot values. "
    "(b) Closed (Complete Genome or Chromosome-level) NCBI assemblies carrying a PL25_M23 BLAST homologue, "
    "broken down by genus and replicon type.",
    ha="center", va="bottom", fontsize=8, color="#555", wrap=True)

plt.tight_layout(rect=[0, 0.04, 1, 1])
for ext in ("pdf","svg","png"):
    fig.savefig(f"{OUTBASE}.{ext}", dpi=300, bbox_inches="tight")
print(f"✓ Saved {OUTBASE}.{{pdf,svg,png}}")
