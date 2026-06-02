# Active-site comparison of PL25_M23 vs four reference M23 enzymes.
# All structures superposed on PL25_M23's M23 catalytic domain (residues 205-330),
# then zoomed to the Zn²⁺-coordinating HxxxD + HxH residues for residue-level comparison.

import os
OUT  = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"
PL25 = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/08_colabfold/PL25_M23_mature/PL25_M23_mature_PL25_M23_mature_form__signal_peptide_removed__residues_38-370__unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"
LYSO = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/4QPB.pdb"
LYTM = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/4ZYB.pdb"
ALE1 = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/ALE1_AF_O05156.pdb"
CWLP = None
import glob
for p in glob.glob("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/14_cwlp_structure/*rank_001*.pdb"):
    CWLP = p; break

# Canonical Zn-binding residues identified from sequence motif scan
SITES = {
    "PL25":  {"obj":"PL25",  "src":PL25, "m23":"205-330", "h1":227, "d1":231, "h2":309, "h3":311, "col":"red"},
    "LYSO":  {"obj":"LYSO",  "src":LYSO, "m23":"all",     "h1":279, "d1":283, "h2":360, "h3":362, "col":"forest"},
    "LYTM":  {"obj":"LYTM",  "src":LYTM, "m23":"all",     "h1":210, "d1":214, "h2":291, "h3":293, "col":"orange"},
    "ALE1":  {"obj":"ALE1",  "src":ALE1, "m23":"140-240", "h1":150, "d1":154, "h2":231, "h3":233, "col":"blue"},
    "CWLP":  {"obj":"CWLP",  "src":CWLP, "m23":"all",     "h1":208, "d1":212, "h2":288, "h3":290, "col":"purple"},
}

cmd.bg_color("white")
cmd.set("antialias", 2)
cmd.set("ray_shadows", 1)
cmd.set("ray_opaque_background", 1)
cmd.set("cartoon_transparency", 0.55)
cmd.set("label_size", 14)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")

# Load and superpose all on PL25 M23 domain
for k, s in SITES.items():
    cmd.load(s["src"], s["obj"])
    cmd.hide("everything", s["obj"])
    cmd.show("cartoon", s["obj"])
    cmd.color(s["col"], s["obj"])

# Superpose each onto PL25 M23 domain
for k, s in SITES.items():
    if k == "PL25": continue
    sel_q = f"{s['obj']} and resi {s['m23']} and name CA+N+C+O" if s["m23"] != "all" else f"{s['obj']} and name CA+N+C+O"
    sel_t = "PL25 and resi 205-330 and name CA+N+C+O"
    try:
        r = cmd.super(s['obj'], "PL25 and resi 205-330", cycles=5, cutoff=2.0)
        print(f"  {k} M23-domain super onto PL25_M23: RMSD={r[0]:.3f} Å, atoms={r[1]}")
    except Exception as e:
        print(f"  ! {k} super failed: {e}")

# Render: show only the Zn-binding side chains as sticks; hide cartoon for clarity
for k, s in SITES.items():
    sel = f"{s['obj']} and resi {s['h1']}+{s['d1']}+{s['h2']}+{s['h3']}"
    cmd.show("sticks", sel)
    cmd.color(s["col"], sel)
    cmd.label(f"{sel} and name CA", "'%s%s' % (resn,resi)")

# Orient on PL25 active site
cmd.orient("PL25 and resi 227+231+309+311")
cmd.zoom("PL25 and resi 227+231+309+311", 4)
cmd.viewport(1800, 1400)
cmd.ray(1800, 1400)
cmd.png(os.path.join(OUT, "09_active_site_overlay.png"), dpi=200)

# Save also as separate per-protein panels (small multiples)
cmd.set("cartoon_transparency", 0.0)
for k, s in SITES.items():
    cmd.disable("all")
    cmd.enable(s["obj"])
    cmd.orient(f"{s['obj']} and resi {s['h1']}+{s['d1']}+{s['h2']}+{s['h3']}")
    cmd.zoom(f"{s['obj']} and resi {s['h1']}+{s['d1']}+{s['h2']}+{s['h3']}", 5)
    cmd.ray(1000, 800)
    cmd.png(os.path.join(OUT, f"09_active_site_{k}.png"), dpi=200)

print("Saved figures/09_active_site_overlay.png + per-protein panels")
