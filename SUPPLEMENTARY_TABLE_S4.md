# Supplementary Table S4

**Independent computational tests of conjugation-system completeness and origin-of-transfer (oriT) presence on PL25 contig_11.**

| Test | Tool / version | Reference profiles searched | Threshold criterion | Result on contig_11 |
|---|---|---|---|---|
| **Plasmid replicon typing** | MOB-suite v3.1.9 | rep_clusters database (Robertson & Nash 2018) | Closest cluster within Mash distance ≤ 0.05 | **No canonical replicon detected** (closest hit: *Spiroplasma citri* CP042474 at Mash distance 0.296, well beyond the typing threshold) |
| **Relaxase typing** | MOB-suite v3.1.9 | MOB family HMMs (MOBP / F / V / Q / C / H) | E ≤ 10⁻³⁰ to any canonical MOB profile | **None of the six canonical MOB families detected** |
| **Mating-pair formation (MPF) typing** | MOB-suite v3.1.9 | MPF_T / F / G / I / B / C / FATA profiles | Standard MPF cluster assignment | **No canonical MPF type assigned** |
| **oriT typing** | MOB-suite v3.1.9 | oriT consensus library | Match to known oriT sequence | **No canonical oriT detected** |
| **Conjugation system completeness** | MacSyFinder v2.1.6 + CONJScan v2.1.0 (Cury et al. 2020) | 17 plasmid models: MOB, T4SS_typeB/C/F/FA/FATA/G/I/T, and decay (dCONJ_*) variants | Mandatory-gene quorum reached per model | **0 complete systems detected** across all 17 models. Partial / sub-quorum hits only: PL25_00148 → T4SS coupling protein (T4SS_t4cp2, HMM score 77); PL25_00150 → VirB11-like ATPase (T4SS_T_virB11, HMM score 79); PL25_00085 → FA_orf13 (weak, score 15). **No MOB-relaxase hit** at any threshold |
| **Canonical oriT motif scan** | Local Python scan (this study) | IncP RP4 (CTGCCTGGTGCT), IncQ RSF1010, IncW R388, IncF consensus, pIP501 Gram⁺ (TAGTGCGCC), MOBP1, ICEBs1 nic-site core | Exact motif match within 2 kb of relaxase candidate + flanking long inverted repeat | **Negative** for all canonical motifs adjacent to the relaxase candidate (PL25_00076). Two chance-level ICEBs1 13-mer matches detected ≥ 24 kb away (positions 1,171 and 30,570 of 74,025 bp) with no inverted-repeat partner. No long (≥ 18 bp) inverted repeat near the relaxase candidate |
| **Origin-of-transfer (oriT) — comprehensive search** | oriTfinder web server (Li et al. 2018) | oriTDB v2 (verified oriT sequences from ICEs, IMEs, conjugative plasmids) | Default thresholds (relaxase + IR + adjacent T4SS) | **To be confirmed** — input file staged at `results/18_oriTfinder_INPUT_contig11.fna`; submit at <https://bioinfo-mml.sjtu.edu.cn/oriTfinder/> (web-only; JavaScript-rendered form prevents programmatic submission) |

## Interpretation

Three independent computational tests (MOB-suite, CONJScan, local oriT-motif scan) converge on the same result: **contig_11 carries divergent HMM-detectable orthologs of individual conjugation-system components (PL25_00076 replication-relaxation protein, PL25_00082 VirB4-like ATPase, PL25_00094 AAA helicase, PL25_00148 T4SS coupling protein) but does not pass canonical thresholds for a complete conjugation system or an origin-of-transfer.** We therefore describe contig_11 throughout the manuscript as a **predicted plasmid-like mobile element carrying conjugation-machinery markers**, rather than as a confirmed conjugative or mobilisable plasmid.

Two non-exclusive interpretations are consistent with these observations:

1. The conjugation system of contig_11 is sufficiently divergent from canonical Gram-positive and Gram-negative reference profiles (built largely from Firmicute/Proteobacterial mesophiles) that the HMMs underlying MOB-suite, CONJScan, and the oriT consensus library fail to detect it at standard thresholds. The strong synteny we observe between contig_11 and three independent NCBI WGS contigs from other halotolerant Bacillaceae genera (Figure 6) is consistent with this interpretation — the lineage exists and is conserved, but it falls in a part of conjugation-system sequence space that current reference databases do not represent.

2. contig_11 represents a non-mobilisable plasmid-like replicon that has retained marker-gene orthologs as evolutionary relics without an intact transfer apparatus.

Distinguishing these requires (i) long-read confirmation of plasmid circularity, (ii) experimental verification of conjugation (mating-out assays into a recipient strain), and (iii) deposition of contig_11 sequence in oriTDB to support future profile updates. All three are flagged as out-of-scope for the present manuscript.

## Methods

- **MOB-suite v3.1.9** was run with `mob_typer --infile contig_11.fna --out_file mob_typer.txt --keep_tmp false` against the v3 reference database (run date 2026-05-17).
- **MacSyFinder v2.1.6 + CONJScan v2.1.0** was run with `--models CONJScan/Plasmids all --db-type ordered_replicon --replicon-topology circular --hmmer hmmsearch -w 4` on the Prokka-predicted protein set of contig_11 (99 CDS).
- **Local oriT motif scan** used in-house Python (Biopython v1.83) to test seven canonical oriT consensus motifs (IncP RP4, IncQ RSF1010, IncW R388, pIP501 Gram⁺, IncF, MOBP1, ICEBs1 nic-site core) on both strands of the 74,025 bp contig_11 sequence, plus a scan for ≥ 18 bp exact inverted repeats within a 6 kb window centred on the predicted relaxase (PL25_00076).
- **oriTfinder** submission was prepared but not executed in this run (web-only tool; manual submission required).
