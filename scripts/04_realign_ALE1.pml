# Re-align ALE-1 onto PL25_M23 using the M23 catalytic domain only.
# The original 10.23 Å RMSD came from superposing onto PDB 1R77, which contains
# ONLY the SH3b cell-wall-targeting domain (residues 271–362) — an unrelated fold.
# Here we use the full-length AlphaFold2 model of ALE-1 (UniProt O05156, v6)
# and restrict the alignment to the M23 catalytic core (residues ~140–240).

PL25  = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/08_colabfold/PL25_M23_mature/PL25_M23_mature_PL25_M23_mature_form__signal_peptide_removed__residues_38-370__unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"
ALE1  = "/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper/results/13_M23_refs/pdb/ALE1_AF_O05156.pdb"

cmd.load(PL25, "PL25")
cmd.load(ALE1, "ALE1")

# PL25_M23 catalytic domain: residues 205-330 (Pfam M23 boundary)
# ALE-1 M23 catalytic domain: residues 140-240 (contains canonical H148/D152 + H229/H231)
cmd.create("PL25_M23",  "PL25 and resi 205-330")
cmd.create("ALE1_M23",  "ALE1 and resi 140-240")

# Sequence-independent superposition with outlier rejection (matches the rest of S2)
r = cmd.super("ALE1_M23", "PL25_M23", cycles=5, cutoff=2.0)
print("=" * 60)
print("ALE-1 M23-domain re-alignment")
print("=" * 60)
print(f"  RMSD (after refinement)        : {r[0]:.3f} Å")
print(f"  Atoms aligned                  : {r[1]}")
print(f"  Refinement cycles              : {r[5]}")
print(f"  RMSD (before refinement)       : {r[3]:.3f} Å")
print(f"  Atoms initial                  : {r[4]}")
print("=" * 60)
