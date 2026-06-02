"""3-way structural superposition figure of PL25_M23 with CwlP (true SLT+M23
reference) and CwlT (architectural contrast — different family: GH25 + NlpC/P60).

Saves: figures/05cc_M23_3way_alignment.{pdf,svg,png}

Run:   /Applications/PyMOL.app/Contents/bin/pymol -cq scripts/05cc_M23_3way_alignment.py
"""
from pymol import cmd
import os, glob

FIG  = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"
PL25 = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/08_colabfold/PL25_M23_mature/PL25_M23_mature_PL25_M23_mature_form__signal_peptide_removed__residues_38-370__unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"
CWLT = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/CwlT_AF_P96645.pdb"   # real CwlT
CWLP = glob.glob("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/14_cwlp_structure/*rank_001*.pdb")[0]

PL25_M23_RES   = "205-330"
PL25_NTERM_RES = "33-185"
# CwlP-clipped numbering (residues 1-314, originals 1373-1686)
CWLP_M23_RES   = "188-314"   # ~last ~125 res = M23 core
CWLP_SLT_RES   = "1-180"
# CwlT (329 aa): GH25 ~40-180 ; NlpC/P60 ~200-329
CWLT_GH25_RES  = "40-180"
CWLT_NLPC_RES  = "200-329"

# Common visualisation settings
def reset():
    cmd.delete("all")
    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_fancy_helices", 1)

def load_all():
    cmd.load(PL25, "PL25")
    cmd.load(CWLT, "CwlT")
    cmd.load(CWLP, "CwlP")
    cmd.hide("everything")
    cmd.show("cartoon", "all")
    cmd.color("red",     "PL25")
    cmd.color("orange",  "CwlT")
    cmd.color("purple",  "CwlP")
    cmd.set("cartoon_transparency", 0.0)

# ---------- Panel a : M23 / C-terminal peptidase domain ONLY ----------
reset(); load_all()
cmd.create("PL25_M23", f"PL25 and resi {PL25_M23_RES}")
cmd.create("CwlT_Cterm", f"CwlT and resi {CWLT_NLPC_RES}")
cmd.create("CwlP_M23", f"CwlP and resi {CWLP_M23_RES}")
cmd.hide("everything"); cmd.show("cartoon", "PL25_M23 or CwlT_Cterm or CwlP_M23")
cmd.color("red", "PL25_M23"); cmd.color("orange","CwlT_Cterm"); cmd.color("purple","CwlP_M23")
r_p = cmd.super("CwlP_M23", "PL25_M23", cycles=5, cutoff=2.0)
r_t = cmd.super("CwlT_Cterm", "PL25_M23", cycles=5, cutoff=2.0)
print(f"  Panel a  PL25-vs-CwlP_M23   : RMSD={r_p[0]:.2f} A ({r_p[1]} atoms)")
print(f"  Panel a  PL25-vs-CwlT_Cterm : RMSD={r_t[0]:.2f} A ({r_t[1]} atoms)")
cmd.orient("PL25_M23"); cmd.zoom("PL25_M23", 4)
cmd.viewport(1400, 1200); cmd.ray(1400, 1200)
cmd.png(f"{FIG}/_05cc_panelA_M23only.png", dpi=200)

# ---------- Panel b : N-terminal glycan-hydrolase domain ONLY ----------
reset(); load_all()
cmd.create("PL25_N",  f"PL25 and resi {PL25_NTERM_RES}")
cmd.create("CwlT_GH25", f"CwlT and resi {CWLT_GH25_RES}")
cmd.create("CwlP_SLT", f"CwlP and resi {CWLP_SLT_RES}")
cmd.hide("everything"); cmd.show("cartoon", "PL25_N or CwlT_GH25 or CwlP_SLT")
cmd.color("red", "PL25_N"); cmd.color("orange","CwlT_GH25"); cmd.color("purple","CwlP_SLT")
r_pn = cmd.super("CwlP_SLT", "PL25_N", cycles=5, cutoff=2.0)
r_tn = cmd.super("CwlT_GH25", "PL25_N", cycles=5, cutoff=2.0)
print(f"  Panel b  PL25_N-vs-CwlP_SLT  : RMSD={r_pn[0]:.2f} A ({r_pn[1]} atoms)")
print(f"  Panel b  PL25_N-vs-CwlT_GH25 : RMSD={r_tn[0]:.2f} A ({r_tn[1]} atoms)")
cmd.orient("PL25_N"); cmd.zoom("PL25_N", 4)
cmd.viewport(1400, 1200); cmd.ray(1400, 1200)
cmd.png(f"{FIG}/_05cc_panelB_Nterm.png", dpi=200)

# ---------- Panel c : Full catalytic core ----------
reset(); load_all()
r_pfull = cmd.super("CwlP", "PL25", cycles=5, cutoff=2.0)
r_tfull = cmd.super("CwlT", "PL25", cycles=5, cutoff=2.0)
print(f"  Panel c  PL25-vs-CwlP FULL : RMSD={r_pfull[0]:.2f} A ({r_pfull[1]} atoms)")
print(f"  Panel c  PL25-vs-CwlT FULL : RMSD={r_tfull[0]:.2f} A ({r_tfull[1]} atoms)")
cmd.orient("PL25"); cmd.zoom("PL25", 3)
cmd.viewport(1400, 1200); cmd.ray(1400, 1200)
cmd.png(f"{FIG}/_05cc_panelC_full.png", dpi=200)

print("\nPanel renders complete. Now run 05cc_combine.py to assemble the final figure.")
