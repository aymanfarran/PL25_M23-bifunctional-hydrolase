# Supplementary Table S3

**Conservation of the canonical M23 Zn²⁺-binding HxxxD + HxH motifs across PL25_M23 and four reference cell-wall hydrolases.**

Residues were identified by sequence motif scan of each parent structure / model and verified by visual inspection of the M23 catalytic site after pairwise superposition onto PL25_M23 (PyMOL v3.1.6.1, `cmd.super`, five outlier-rejection cycles at 2.0 Å). M23-domain RMSDs (Cα) versus PL25_M23 residues 205–330 are listed for context; all five enzymes superpose into a narrow 0.42–0.68 Å band, confirming a structurally conserved catalytic core.

| Protein (parent structure) | Source / location | HxxxD: H₁ | HxxxD: D | HxH: H₂ | HxH: H₃ | M23 RMSD vs PL25 (Å) | Atoms |
|---|---|---:|---:|---:|---:|---:|---:|
| **PL25_M23** (this study, AF2) | *V. salarius* PL25 plasmid contig_11 | **H227** | **D231** | **H309** | **H311** | — | — |
| Lysostaphin (PDB 4QPB) | *Staphylococcus simulans*, secreted enzybiotic | H279 | D283 | H360 | H362 | 0.50 | 486 |
| LytM (PDB 4ZYB) | *Staphylococcus aureus*, autolysin | H210 | D214 | H291 | H293 | 0.60 | 546 |
| ALE-1 (AF2 of UniProt O05156)¹ | *Staphylococcus capitis*, secreted lysin | H150 | D154 | H231 | H233 | 0.42 | 398 |
| CwlP-clipped (AF2 of UniProt O31976, residues 1373–1686)² | *Bacillus subtilis*, SP-β prophage; reported Zn²⁺-independent³ | H1580 (H208 of clip) | D1584 (D212) | H1660 (H288) | H1662 (H290) | 0.68 | 540 |

## Notes

1. **ALE-1**: PDB 1R77 deposited for this protein covers only the C-terminal SH3b cell-wall-targeting domain (residues 271–362) and does not include the M23 catalytic domain. The AlphaFold2 monomer (AF-O05156-F1, v6) is used here as it provides full-length coverage including the M23 catalytic core (residues ~140–240); this corrects an earlier superposition error (10.23 Å) that arose from comparing PL25's M23 fold against ALE-1's SH3b fold.
2. **CwlP-clipped**: full-length CwlP is a 2,285-aa SP-β prophage tail-associated protein; the SLT + M23 catalytic core (residues 1373–1686) was modelled with ColabFold v1.5.5 (AlphaFold2 ptm). Residue numbering above gives both the original full-length position and (in parentheses) the position in the 314-aa clipped construct.
3. **CwlP Zn²⁺-independence**: Sudiarta et al. (2010, *J. Biol. Chem.* 285, 41232–41243) reported that CwlP's M23 domain does not require Zn²⁺ for activity. Importantly, our sequence motif scan shows that CwlP **retains** both canonical motifs at the primary-sequence level (H1580–D1584 and H1660/H1662); the reported Zn²⁺-independence therefore cannot be inferred from the motif alone and must reflect more subtle features (e.g. coordination geometry, second-shell residues, substrate-induced conformational requirements). Any structural claim of Zn²⁺-dependence for PL25_M23 made in the main text is therefore based on **conservation of the canonical motif** rather than on a CwlP-style motif loss; experimental verification (EDTA inhibition, Zn²⁺ titration) is required to confirm Zn²⁺-dependence biochemically.
4. **CwlT**: PDB 4FDY, used in earlier analyses as a "CwlT (*B. subtilis* ICEBs1)" reference, actually maps to UniProt Q932I6 (a *Staphylococcus aureus* NlpC/P60 family lipoprotein) and is **not** *B. subtilis* CwlT (UniProt O32008). CwlT is therefore omitted from this active-site table pending re-identification of an appropriate experimental or AF2 structure. All structural comparisons involving "CwlT" in Figures 5 and Supplementary Table S2 should be re-checked against this potential mis-assignment before submission.

## Methods

- AlphaFold2 models obtained from the EBI AlphaFold Database (v6).
- Motif scan: in-house Python script searching the primary sequence of each parent structure for the regular-expression patterns `H[A-Z]{3}D` (HxxxD) and `H[A-Z]H` (HxH) within the predicted M23 catalytic domain.
- Pairwise structural superposition: PyMOL v3.1.6.1, `cmd.super` with five cycles of outlier rejection at 2.0 Å; reference selection PL25_M23 residues 205–330.
- All side chains and Cα atoms used in the alignment table are visualised in main-text Figure 9 (active-site overlay).
