#!/usr/bin/env python3
"""M23 phylogeny v2 — PL25 + 6 refs + 14 phage M23 + 80 Bacillaceae homologues.
Colour-coded by source: PL25 / Reference / Phage / Host genus."""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

ROOT = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
TREE = ROOT/"results/07_M23_phylogeny/M23_iqtree_v2.treefile"
OUT  = ROOT/"figures"

# Color scheme by category
CAT_COLORS = {
    "PL25":              "#D7191C",
    "Reference":         "#9C27B0",
    "Phage_M23":         "#FF7F00",      # orange — phage M23 lysins (NEW category)
    "Virgibacillus":     "#FF7F0E",
    "Bacillus":          "#1F77B4",
    "Lentibacillus":     "#2CA02C",
    "Cytobacillus":      "#E377C2",
    "Ornithinibacillus": "#FFBB78",
    "Halophile_other":   "#B0BEC5",
    "Caldifermentibacillus":"#FF9800",
    "Caldibacillus":     "#FFA726",
    "Domibacillus":      "#7E57C2",
    "Niallia":           "#26A69A",
    "Other":             "#999999",
}

tree = Phylo.read(str(TREE), "newick")
try: tree.root_at_midpoint()
except: pass
n_terms = len(tree.get_terminals())
print(f"Tree leaves: {n_terms}")

# Build taxonomy / category lookup from seed set FASTA
tax_map = {}
with open(ROOT/"results/07_M23_phylogeny/M23_seedset_v2.faa") as f:
    for line in f:
        if not line.startswith(">"): continue
        head = line[1:].strip()
        ide = head.split()[0]
        if ide.startswith("PL25"):
            tax_map[ide] = "PL25"
        elif ide.startswith("REF_"):
            tax_map[ide] = "Reference"
        elif ide.startswith("PHAGE_"):
            tax_map[ide] = "Phage_M23"
        elif "[" in head and "]" in head:
            organism = head[head.rfind("[")+1:head.rfind("]")]
            genus = organism.split()[0]
            tax_map[ide] = genus

# Re-categorise
def cat_for(name):
    if name in tax_map:
        c = tax_map[name]
        if c in CAT_COLORS: return c
        halo = {"Aquisalibacillus","Oceanobacillus","Piscibacillus","Halobacillus",
                "Gracilibacillus","Ureibacillus","Margalitia","Solibacillus","Peribacillus",
                "Siminovitchia","Lysinibacillus","Carnobacterium","Streptomyces","Paenibacillus"}
        if c in halo: return "Halophile_other"
    return "Other"

# Coords
leaves = tree.get_terminals()
y = {id(lf): i for i, lf in enumerate(leaves)}
def assign_y(c):
    if c.is_terminal(): return y[id(c)]
    ys = [assign_y(ch) for ch in c.clades]
    y[id(c)] = (min(ys)+max(ys))/2
    return y[id(c)]
assign_y(tree.root)
x = {}
def assign_x(c, d=0.0):
    x[id(c)] = d
    for ch in c.clades:
        assign_x(ch, d + (ch.branch_length or 0.0))
assign_x(tree.root)
xmax = max(x.values())

# Figure
fig_h = max(15, 0.20*n_terms + 3)
fig, ax = plt.subplots(figsize=(14, fig_h))

# Branches
def draw(c):
    cy, cx = y[id(c)], x[id(c)]
    for ch in c.clades:
        chy, chx = y[id(ch)], x[id(ch)]
        ax.plot([cx, cx], [cy, chy], color="#444", lw=0.8, zorder=2)
        ax.plot([cx, chx], [chy, chy], color="#444", lw=0.8, zorder=2)
        draw(ch)
draw(tree.root)

# UFBoot labels
for c in tree.get_nonterminals():
    if c is tree.root: continue
    bs = c.confidence
    if bs is None:
        try: bs = float(c.name)
        except: bs = None
    if bs is None or bs < 70: continue
    cx, cy = x[id(c)], y[id(c)]
    if bs >= 95: bc, fw = "#1A7A1A", "bold"
    elif bs >= 85: bc, fw = "#888800", "bold"
    else: bc, fw = "#666666", "normal"
    ax.text(cx-xmax*0.003, cy, f"{int(bs)}",
            ha="right", va="center", fontsize=5,
            color=bc, fontweight=fw, zorder=4,
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.7))

# Leaves
for lf in leaves:
    ly, lx = y[id(lf)], x[id(lf)]
    cat = cat_for(lf.name)
    col = CAT_COLORS.get(cat, "#999999")
    name = lf.name
    if name.startswith("REF_"):
        label = name.replace("REF_","").replace("_"," ")[:40]
        fs, fw, lcol = 7.5, "bold", "#9C27B0"
        marker, ms = "s", 50
    elif name.startswith("PL25"):
        label = "★  PL25 M23  (this study)"
        fs, fw, lcol = 9.5, "bold", "#D7191C"
        marker, ms = "*", 110
    elif name.startswith("PHAGE_"):
        label = name.replace("PHAGE_","").replace("_"," ")[:45]
        fs, fw, lcol = 7.5, "bold", "#FF6F00"
        marker, ms = "D", 50  # diamond
    else:
        g = tax_map.get(name, "")
        label = f"{name}  [{g}]"[:42]
        fs, fw, lcol = 6.0, "normal", "#222"
        marker, ms = "o", 16
    ax.text(lx+xmax*0.005, ly, label,
            va="center", ha="left", fontsize=fs,
            fontweight=fw, color=lcol, zorder=4)
    ax.scatter([lx], [ly], s=ms, marker=marker, color=col,
               edgecolor="black", linewidths=0.5, zorder=5)

ax.set_xlim(-xmax*0.04, xmax*1.40)
ax.set_ylim(-1, n_terms)
ax.invert_yaxis()
ax.set_xlabel("Substitutions per site")
ax.set_yticks([])
for s in ("top","right","left"): ax.spines[s].set_visible(False)
ax.grid(axis="x", alpha=0.20)

ax.set_title("M23 metallopeptidase ML phylogeny — PL25 conjugative plasmid lysin in context\n"
             "IQ-TREE3 · Q.PFAM+G4 · 1,000 UFBoot · 101 sequences · midpoint-rooted",
             fontsize=12, fontweight="bold", loc="left", pad=10)

# Footnote
fig.text(0.5, 0.011,
    "★ red = PL25 M23 (this study, plasmid-encoded)  •  ▢ purple = characterised reference enzymes (lysostaphin, LytM, ALE-1, CwlT, EnpA, ZoocinA)  •  "
    "◆ orange = phage-encoded M23 lysins (HydH5, LysK, PlyGRCS, Twort, A118, B4, etc.)  •  ● = Bacillaceae homologues from NCBI nr.  "
    "UFBoot ≥95 dark green, 85-94 olive, 70-84 grey.",
    ha="center", va="bottom", fontsize=7.5, style="italic", color="#444", wrap=True)

# Legend
seen = set()
ordered = ["PL25","Reference","Phage_M23","Virgibacillus","Lentibacillus","Caldifermentibacillus",
           "Cytobacillus","Ornithinibacillus","Domibacillus","Niallia","Halophile_other","Other"]
items = []
for name in ordered:
    if name in CAT_COLORS:
        items.append(mpatches.Patch(color=CAT_COLORS[name], label=name.replace("_"," ")))
items.append(mpatches.Patch(color="#1A7A1A", label="UFBoot ≥ 95"))
items.append(mpatches.Patch(color="#888800", label="UFBoot 85-94"))
items.append(mpatches.Patch(color="#666666", label="UFBoot 70-84"))
ax.legend(handles=items, loc="lower right", bbox_to_anchor=(1.0, 0.0),
          ncol=2, fontsize=8, title="Category / support",
          title_fontsize=9)

plt.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(OUT/"02b_M23_phylogeny_v2.pdf", bbox_inches="tight")
fig.savefig(OUT/"02b_M23_phylogeny_v2.svg", bbox_inches="tight")
fig.savefig(OUT/"02b_M23_phylogeny_v2.png", dpi=200, bbox_inches="tight")
print(f"✓ Saved 02b_M23_phylogeny_v2.{{pdf,svg,png}}")
