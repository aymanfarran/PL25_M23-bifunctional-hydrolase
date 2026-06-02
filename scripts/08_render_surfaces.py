"""Re-render electrostatic-style surfaces for the 5 proteins under a SINGLE
common color rule, so the panels are directly comparable.

Rule (applied identically to every protein):
  acidic side chains (D, E)            → red    (charge ≈ −1 at pH 7)
  basic side chains  (K, R)            → blue   (charge ≈ +1 at pH 7)
  His (partial at pH 7)                → light pink
  all other residues                   → white

Surface is then coloured by the underlying atom colors, giving a per-atom
charge representation that does not depend on a protein-specific numeric
electrostatic-potential range (the previous panels each had a different
±kT/e scale, which weakens direct visual comparison).

Two views per protein: front and 180° rotation.
Rendered with PyMOL Open-Source (no evaluation watermark).

Run: /usr/local/Caskroom/miniforge/base/bin/pymol -cq scripts/08_render_surfaces.py
"""
from pymol import cmd
import os, glob

OUT  = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/figures"
PL25 = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/08_colabfold/PL25_M23_mature/PL25_M23_mature_PL25_M23_mature_form__signal_peptide_removed__residues_38-370__unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"
LYSO = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/4QPB.pdb"
LYTM = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/4ZYB.pdb"
CWLT = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/CwlT_AF_P96645.pdb"
CWLP = glob.glob("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/14_cwlp_structure/*rank_001*.pdb")[0]

PROTEINS = [
    ("PL25_M23",    PL25),
    ("CwlP",        CWLP),
    ("CwlT",        CWLT),
    ("Lysostaphin", LYSO),
    ("LytM",        LYTM),
]

cmd.bg_color("white")
cmd.set("ray_shadows", 0)
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("surface_quality", 1)

# Common color rule
COL_ACID  = "tv_red"
COL_BASE  = "tv_blue"
COL_HIS   = "lightpink"
COL_NEUT  = "white"

for name, src in PROTEINS:
    cmd.delete("all")
    cmd.load(src, name)
    # Choose chain A only (avoid surfaces of crystal dimers)
    cmd.remove(f"{name} and not chain A")
    # If chain A is empty (single-chain AF2), no-op
    if cmd.count_atoms(name) == 0:
        cmd.load(src, name)

    cmd.hide("everything")
    cmd.show("surface", name)
    # Default everything white
    cmd.color(COL_NEUT, name)
    # Acidic side-chain atoms
    cmd.color(COL_ACID, f"{name} and resn ASP+GLU and not name N+CA+C+O")
    # Basic side-chain atoms
    cmd.color(COL_BASE, f"{name} and resn LYS+ARG and not name N+CA+C+O")
    # His
    cmd.color(COL_HIS,  f"{name} and resn HIS and not name N+CA+C+O")

    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.0)

    # Front view
    cmd.orient(name)
    cmd.zoom(name, 2.5)
    cmd.viewport(900, 900)
    cmd.ray(900, 900)
    cmd.png(os.path.join(OUT, f"_es2_{name}_front.png"), dpi=200)

    # 180° rotation
    cmd.rotate("y", 180)
    cmd.zoom(name, 2.5)
    cmd.ray(900, 900)
    cmd.png(os.path.join(OUT, f"_es2_{name}_back.png"), dpi=200)

print("Surface panels rendered (no watermark, common color rule).")
