# Supplementary Table S1

**Assembly metrics, Flye-reported circularity, and geNomad classification of all PL25 contigs.**

The 4.75 Mb *Virgibacillus salarius* PL25 genome was assembled with Flye v2.9.3-b1797 from Oxford Nanopore PromethION reads (`--nano-raw`). Flye-reported per-contig metrics are taken from `assembly_info.txt`: `circ` indicates the contig was reported as a closed circle by the assembler (independent confirmation by junction-spanning-read mapping was not performed); `cov.` is the mean mapped coverage; `mult.` is Flye's multiplicity estimate. geNomad v1.7.0 was run end-to-end with relaxed marker-enrichment thresholds (`--min-virus-marker-enrichment 0 --min-plasmid-marker-enrichment 0 --min-virus-hallmarks 0 --min-plasmid-hallmarks 0`) to improve detection of potentially divergent halotolerant mobile elements; default ML score thresholds (min_score 0.7, max_fdr 0.1) were retained. GC1/GC2/GC3 values were computed in-house from Prokka-predicted CDS; Δ-values are differences relative to the chromosomal baseline (contig_4). The PL25_M23 hydrolase (PL25_00078) is encoded on contig_11.

Working-FASTA contig nomenclature was preserved from `PL25.fasta` and differs from the original Flye output for two contigs: the 74,025-bp closed circular plasmid candidate is contig_12 in Flye but contig_11 in the working FASTA, and the 44,498-bp closed circular element is contig_11 in Flye but contig_12 in the working FASTA. Throughout the manuscript and this table, **contig_11** refers to the 74,025-bp plasmid candidate. The full mapping is given in Supplementary Methods S1.

## Working FASTA contigs

| Manuscript contig | Flye contig | Length (bp) | n CDS | Flye `circ` | Flye cov. (×) | Flye mult. | GC % | GC1 % | GC2 % | GC3 % | ΔGC vs chr. | ΔGC3 vs chr. | geNomad plasmid score | geNomad virus score | Final assignment |
|---|---|---:|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **contig_4** | contig_4 | 4,349,560 | 4,098 | N | 476 | 1 | 36.81 | 48.95 | 33.96 | 28.71 | — | — | 0.097 | — | Main chromosomal contig (not closed by Flye) |
| contig_2 | contig_2 | 249,289 | — | N | 849 | 2 | 36.29 | — | — | — | −0.52 | — | — | — | Unannotated chromosomal fragment (no Prokka CDS calls; high coverage and elevated multiplicity are consistent with a chromosomal repeat) |
| **contig_11** | **contig_12** | **74,025** | **99** | **Y** | **889** | **2** | **31.22** | **42.62** | **30.54** | **23.19** | **−5.59** | **−5.52** | **0.989** | — | **Closed circular plasmid candidate (carries PL25_M23)** |
| contig_12 | contig_11 | 44,498 | 86 | Y | 22 | 1 (alt 2) | 36.31 | 45.86 | 32.47 | 30.85 | −0.50 | +2.14 | — | 0.979 | Closed circular element, repeat-flagged; geNomad virus-classified (likely a free / extrachromosomal prophage form) |
| contig_10 | contig_10 | 32,390 | 70 | Y | 86 | 1 | 29.70 | 41.04 | 29.63 | 20.59 | −7.11 | −8.12 | 0.967 | — | Closed circular replicon |
| contig_9 | contig_9 | 3,944 | 5 | N | 4 | 1 | 35.45 | 45.71 | 33.13 | 26.99 | −1.36 | −1.72 | 0.818 | — | Short linear fragment (geNomad plasmid score above threshold but unclosed) |
| contig_6 | contig_6 | 514 | 0 | N | 3 | 1 | 40.86 | — | — | — | +4.05 | — | — | — | Short linear fragment (< 600 bp; no Prokka CDS) |
| contig_5 | contig_5 | 462 | 0 | N | 4 | 4 | 41.34 | — | — | — | +4.53 | — | — | — | Short linear fragment (< 600 bp; no Prokka CDS) |

## Contig excluded from the working FASTA

| Flye contig | Length (bp) | Flye `circ` | Flye cov. (×) | Reason for exclusion |
|---|---:|:---:|---:|---|
| contig_1 | 3,815 | N | 4 | Failed the pre-defined inclusion criterion (minimum coverage ≥ 10× and minimum length ≥ 5 kb). Excluded from the working FASTA and from the GenBank submission. |

## Provirus regions detected within the chromosome (contig_4)

Six high-confidence provirus regions were predicted by geNomad inside contig_4 (chromosome). These are integrated and so are not reported as separate contigs in the per-contig table above; they are listed here for completeness.

| Region (coords on contig_4) | Length (bp) | n CDS | virus score | n hallmarks | Taxonomy |
|---|---:|---:|---:|---:|---|
| 735,616–777,289 | 41,674 | 81 | 0.9791 | 13 | Caudoviricetes |
| 506,914–544,733 | 37,820 | 51 | 0.9788 | 12 | Caudoviricetes |
| 1,208,625–1,251,666 | 43,042 | 63 | 0.9696 | 11 | Caudoviricetes |
| 890,218–930,103 | 39,886 | 51 | 0.9690 | 10 | Caudoviricetes |
| 2,883,721–2,896,061 | 12,341 | 15 | 0.9188 | 4 | Caudoviricetes |

## Notes

- **Sequencing and assembly.** The Nanopore PromethION input (496,746 reads; 2.40 Gb; N50 10,372 bp; mean read length 4,830 bp; average Phred Q17.5; Q20 fraction 87.5 %) provided ~510× theoretical coverage of the 4.7 Mb genome. Assembly used Flye v2.9.3-b1797 with `--nano-raw`; no external polishing step was applied. Circularity of individual contigs was inferred from the Flye `assembly_info.txt` (`circ` field) and was not independently verified by junction-spanning-read mapping.
- **CDS calling.** Prokka v1.15.6 (Prodigal v2.6.3) was used with `--kingdom Bacteria --genus Virgibacillus --species salarius --strain PL25 --rfam --addgenes --locustag PL25` to obtain 4,613 CDS across the eight working-FASTA contigs. GC1/GC2/GC3 were calculated by an in-house Python script counting GC at each codon position across all predicted CDS per contig.
- **Chromosome status.** contig_4 was not reported as a closed circle by Flye; it is therefore described as the *main chromosomal contig*, not as a "closed chromosome." Independent chromosomal-closure evidence (e.g. long-read junction support, optical mapping) has not been collected.
- **Interpretation of compositional offsets.** Negative ΔGC and ΔGC3 values larger than the within-chromosome variance are consistent with a distinct evolutionary history and possible horizontal acquisition (Lawrence & Ochman, *PNAS* 1998). Such compositional shifts do not establish recent transfer; they must be interpreted alongside the conjugation-machinery and synteny analyses of Supplementary Table S4 and Figure 5.
