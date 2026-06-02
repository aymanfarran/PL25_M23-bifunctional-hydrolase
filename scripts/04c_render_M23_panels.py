"""Re-render the four 3D structure views of PL25_M23 with FULL-LENGTH residue
labels (H264, D268, H346, H348). PDB uses mature numbering 1–333; offset = +37.

Run:  /Applications/PyMOL.app/Contents/bin/pymol -cq scripts/04c_render_M23_panels.py
"""
from pymol import cmd
import os

OUT = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"
PDB = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/08_colabfold/PL25_M23_mature/PL25_M23_mature_PL25_M23_mature_form__signal_peptide_removed__residues_38-370__unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"

# Domain ranges in MATURE numbering (PDB residue numbers). +37 = full-length.
DOM_SLT35_M   = (33, 185)     # full-length 70–222
DOM_LINKER_M  = (186, 225)    # full-length 223–262
DOM_M23_M     = (226, 323)    # full-length 263–360
ZN_M          = [227, 231, 309, 311]  # full-length 264, 268, 346, 348

cmd.load(PDB, "M23")
cmd.bg_color("white")
cmd.hide("everything")
cmd.show("cartoon", "M23")
cmd.remove("M23 and resi 1-30")  # trim disordered N-terminal residues

cmd.color("forest",     f"M23 and resi {DOM_SLT35_M[0]}-{DOM_SLT35_M[1]}")
cmd.color("paleyellow", f"M23 and resi {DOM_LINKER_M[0]}-{DOM_LINKER_M[1]}")
cmd.color("firebrick",  f"M23 and resi {DOM_M23_M[0]}-{DOM_M23_M[1]}")

cmd.select("active_site", f"resi {'+'.join(map(str, ZN_M))}")
cmd.show("sticks", "active_site")
cmd.color("orange", "active_site")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("antialias", 2)
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 1)
cmd.set("label_size", 18)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")

# Front
cmd.orient("M23"); cmd.zoom("M23", 3)
cmd.viewport(1400, 1200); cmd.ray(1400, 1200)
cmd.png(os.path.join(OUT, "_panel_A_front_v2.png"), dpi=200)

# 90° rotation
cmd.rotate("y", 90); cmd.zoom("M23", 3)
cmd.ray(1400, 1200); cmd.png(os.path.join(OUT, "_panel_B_side_v2.png"), dpi=200)

# 180° back view
cmd.rotate("y", 90); cmd.zoom("M23", 3)
cmd.ray(1400, 1200); cmd.png(os.path.join(OUT, "_panel_C_back_v2.png"), dpi=200)

# Active site close-up with full-length residue labels (mature+37)
cmd.orient("active_site"); cmd.zoom("active_site", 5)
cmd.label("active_site and name CA",
          "'%s%s' % (resn, str(int(resi)+37))")    # ← full-length numbering
cmd.ray(1400, 1200); cmd.png(os.path.join(OUT, "_panel_D_activesite_v2.png"), dpi=200)

print("Re-rendered 4 PyMOL panels with full-length residue labels.")
