"""Active-site overlay + per-protein close-ups for PL25_M23 vs four references.
Run with:  /Applications/PyMOL.app/Contents/bin/pymol -cq scripts/05_active_site_panel.py
"""
from pymol import cmd
import os, glob

OUT  = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"
PL25 = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/08_colabfold/PL25_M23_mature/PL25_M23_mature_PL25_M23_mature_form__signal_peptide_removed__residues_38-370__unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"
LYSO = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/4QPB.pdb"
LYTM = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/4ZYB.pdb"
ALE1 = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/ALE1_AF_O05156.pdb"
CWLP = glob.glob("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/14_cwlp_structure/*rank_001*.pdb")[0]

SITES = [
    {"k":"PL25", "src":PL25, "m23":"205-330", "h1":227, "d1":231, "h2":309, "h3":311, "col":"red",     "label":"PL25_M23 (this study, AF2)"},
    {"k":"LYSO", "src":LYSO, "m23":None,      "h1":279, "d1":283, "h2":360, "h3":362, "col":"forest",  "label":"Lysostaphin (4QPB)"},
    {"k":"LYTM", "src":LYTM, "m23":None,      "h1":210, "d1":214, "h2":291, "h3":293, "col":"orange",  "label":"LytM (4ZYB)"},
    {"k":"ALE1", "src":ALE1, "m23":"140-240", "h1":150, "d1":154, "h2":231, "h3":233, "col":"blue",    "label":"ALE-1 (AF2, O05156)"},
    {"k":"CWLP", "src":CWLP, "m23":None,      "h1":208, "d1":212, "h2":288, "h3":290, "col":"purple",  "label":"CwlP-clipped (AF2)"},
]

cmd.bg_color("white")
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("ray_opaque_background", 1)
cmd.set("cartoon_transparency", 0.6)
cmd.set("label_size", 14)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")
cmd.set("stick_radius", 0.22)

for s in SITES:
    cmd.load(s["src"], s["k"])
    cmd.hide("everything", s["k"])
    cmd.show("cartoon", s["k"])
    cmd.color(s["col"], s["k"])

# Superpose each reference onto PL25 M23 catalytic domain
for s in SITES:
    if s["k"] == "PL25": continue
    sel_q = f"{s['k']} and resi {s['m23']}" if s["m23"] else s["k"]
    r = cmd.super(sel_q, "PL25 and resi 205-330", cycles=5, cutoff=2.0)
    print(f"  {s['k']:6s} super -> PL25_M23: RMSD={r[0]:.3f} Å, atoms={r[1]}")

# Show catalytic side chains as sticks
for s in SITES:
    sel = f"{s['k']} and resi {s['h1']}+{s['d1']}+{s['h2']}+{s['h3']}"
    cmd.show("sticks", sel)
    cmd.color(s["col"], sel)
    cmd.label(f"{sel} and name CA", "'%s%s' % (resn,resi)")

# Render overlay
cmd.orient("PL25 and resi 227+231+309+311")
cmd.zoom("PL25 and resi 227+231+309+311", 5)
cmd.viewport(1800, 1400)
cmd.ray(1800, 1400)
cmd.png(os.path.join(OUT, "09_active_site_overlay.png"), dpi=200)
print("✓ Saved 09_active_site_overlay.png")

# Per-protein close-ups
cmd.set("cartoon_transparency", 0.0)
for s in SITES:
    cmd.disable("all")
    cmd.enable(s["k"])
    sel = f"{s['k']} and resi {s['h1']}+{s['d1']}+{s['h2']}+{s['h3']}"
    cmd.orient(sel)
    cmd.zoom(sel, 6)
    cmd.ray(1000, 800)
    cmd.png(os.path.join(OUT, f"09_active_site_{s['k']}.png"), dpi=200)
    print(f"  saved 09_active_site_{s['k']}.png")

# Combine into a single 5-panel figure with matplotlib
cmd.disable("all")
print("Done — now run scripts/05b_combine_panels.py to assemble the 5-panel figure.")
