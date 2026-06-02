# Supplementary Table S2

**Structural superposition of PL25_M23 against reference cell-wall hydrolases.**
All pairwise superpositions were performed in PyMOL 2.5.5 using the `super` command (sequence-independent alignment with five cycles of outlier rejection at 2.0 Å). The query structure was the AlphaFold2 ptm model of mature PL25_M23 (residues 38–370, rank 1, mean pLDDT = 88.7). Reference structures are experimental crystal structures from the PDB unless otherwise noted; CwlP-clipped is an AlphaFold2 model of the catalytic core (residues 1373–1686) of the 2,285-aa SP-β prophage protein (UniProt O31976). "Atoms aligned" is the number of Cα atom pairs retained after outlier rejection.

## A. M23 catalytic domain (PL25 residues 205–330) versus reference M23 domains

| Reference | Source organism | PDB / model | Function | RMSD (Å) | Atoms aligned |
|---|---|---|---|---:|---:|
| **CwlP M23 (clipped)** | *Bacillus subtilis* (SP-β prophage) | AF2 model (UniProt O31976, residues 1373–1686) | Prophage-encoded M23 (Zn²⁺-independent) | **0.63** | 466 |
| **Lysostaphin** | *Staphylococcus simulans* | 4QPB | Anti-MRSA enzybiotic | **0.49** | 486 |
| **LytM** | *Staphylococcus aureus* | 4ZYB | Autolysin | **0.60** | 550 |
| **EnpA** | *Enterococcus faecalis* (prophage) | 6S6E | Phage-encoded M23 | 0.77 | 36 |
| **ZoocinA** | *Streptococcus equi* subsp. *zooepidemicus* | 5KVP | Anti-streptococcus M23 | 1.03 | 458 |
| **CwlT M23 (clipped)** | *Bacillus subtilis* (ICE *ICEBs1*) | 4FDY (chain A, M23 domain only) | ICE-encoded conjugation lysin | 1.17 | 65 |
| ALE-1 | *Staphylococcus capitis* | 1R77 | Anti-staphylococcus M23 | 10.23 | 318 (poor; different chain orientation) |

## B. N-terminal glycoside-hydrolase domain (PL25 residues 33–185) versus reference glycoside-hydrolase domains

| Reference | Source / domain | Fold | RMSD (Å) | Atoms aligned |
|---|---|---|---:|---:|
| **CwlP SLT (clipped)** | *B. subtilis* SP-β prophage, SLT (PF01464) domain | Lysozyme-like / Slt35-like | **1.31** | 401 |
| CwlT GH25 | *B. subtilis* ICE *ICEBs1*, GH25 muramidase domain | GH25 (CAZy) — different fold | 7.92 | 490 (poor; non-homologous fold) |

## C. Full bifunctional catalytic core (SLT + M23)

| Comparison | Architecture | RMSD (Å) | Notes |
|---|---|---:|---|
| PL25_M23 vs CwlP-clipped | SLT + M23 (same architecture) | reported domain-wise (see A & B) | Inter-domain orientation differs; superposition forced on individual domains |
| PL25_M23 vs CwlT (4FDY) | GH25 + M23 (different N-terminal fold) | not reported as a global RMSD | N-terminal folds non-homologous; only M23 domain comparable |

## Methods notes

- **Software:** PyMOL 2.5.5 (`super` command); cycles = 5; outlier rejection = 2.0 Å; gap-open and gap-extension penalties at PyMOL defaults.
- **PL25_M23 model:** ColabFold v1.5.5 (AlphaFold2 ptm + 3 recycles); top-ranked model selected; mean pLDDT = 88.7; 70.6 % of residues at pLDDT ≥ 90.
- **CwlP-clipped model:** ColabFold v1.5.5 run on UniProt O31976 residues 1373–1686 (covering SLT + M23 domains identified by InterProScan 5.61-93.0); mean pLDDT > 85 across both catalytic domains.
- **PDB references:** 4QPB (lysostaphin), 4ZYB (LytM), 6S6E (EnpA), 5KVP (ZoocinA), 4FDY (CwlT), 1R77 (ALE-1).
- All atoms reported are Cα atoms retained after PyMOL's iterative outlier-rejection cycles.
- "Domain-only" comparisons used the PyMOL selection limited to the indicated residue ranges before invoking `super`.

The values reported in this table are summarised in the main-text Figure 5 (PL25 vs CwlT vs CwlP-clipped three-way superposition) and Figures 5a/5b (six- and three-reference panels, deposited as supplementary figures).
