# Supplementary Table S1

**geNomad 1.7.0 classification and GC composition of all PL25 contigs.**
Plasmid and virus scores are posterior probabilities from the geNomad ML classifier (range 0–1; > 0.7 considered confident). GC3 is the GC content at the third codon position (computed by an in-house script over Prokka-predicted CDS) and is the most sensitive single indicator of recent horizontal acquisition. Δ-values are differences relative to the chromosomal baseline (contig_4). The PL25 M23 hydrolase (PL25_00078) is encoded on contig_11.

| Contig | Length (bp) | n CDS | GC % | GC1 % | GC2 % | GC3 % | ΔGC vs chr. | ΔGC3 vs chr. | geNomad plasmid score | geNomad virus score | Final assignment |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **contig_4** | 4,349,560 | 4,098 | 36.81 | 48.95 | 33.96 | 28.71 | — | — | 0.097 | — | Chromosome (baseline) |
| **contig_11** | 74,025 | 99 | 31.22 | 42.62 | 30.54 | 23.19 | **−5.59** | **−5.52** | **0.989** | — | **Plasmid-like element (M23 carrier)** |
| contig_10 | 32,390 | 70 | 29.70 | 41.04 | 29.63 | 20.59 | −7.11 | −8.12 | 0.967 | — | Plasmid-like element |
| contig_9 | 3,944 | 5 | 35.45 | 45.71 | 33.13 | 26.99 | −1.36 | −1.72 | 0.818 | — | Small plasmid-like element |
| contig_12 | 44,498 | 86 | 36.31 | 45.86 | 32.47 | 30.85 | −0.50 | +2.14 | — | 0.979 | Virus / prophage-like (free) |
| contig_2 | 249,289 | 0 | 36.29 | — | — | — | −0.52 | — | — | — | Unannotated chromosomal fragment |
| contig_5 | 462 | 0 | 41.34 | — | — | — | +4.53 | — | — | — | Short contig (no CDS) |
| contig_6 | 514 | 0 | 40.86 | — | — | — | +4.05 | — | — | — | Short contig (no CDS) |

## Provirus regions detected in the chromosome (contig_4)

Six high-confidence provirus regions were predicted by geNomad within contig_4 (chromosome) but are not reported in the per-contig table above because they are integrated:

| Region (coords on contig_4) | Length (bp) | n CDS | virus score | n hallmarks | Taxonomy |
|---|---:|---:|---:|---:|---|
| 735,616–777,289 | 41,674 | 81 | 0.9791 | 13 | Caudoviricetes |
| 506,914–544,733 | 37,820 | 51 | 0.9788 | 12 | Caudoviricetes |
| 1,208,625–1,251,666 | 43,042 | 63 | 0.9696 | 11 | Caudoviricetes |
| 890,218–930,103 | 39,886 | 51 | 0.9690 | 10 | Caudoviricetes |
| 2,883,721–2,896,061 | 12,341 | 15 | 0.9188 | 4 | Caudoviricetes |

## Notes

- **geNomad v1.7.0** was run in `end-to-end` mode (`--cleanup --splits 4`) with default ML thresholds (min_score 0.7, max_fdr 0.1) on the unmodified Unicycler assembly `PL25.fasta` (MD5 `fe7d19423bf833399f1222bbcd184195`).
- **CDS were predicted with Prokka v1.15.6** (--genus *Virgibacillus*); GC1/GC2/GC3 were calculated by an in-house Python script counting GC at each codon position across all predicted CDS per contig.
- The chromosomal baseline (contig_4) GC = 36.81 % and GC3 = 28.71 % is used as the reference for Δ-values.
- Negative ΔGC3 values larger than the within-chromosome variance are interpreted as evidence of recent horizontal acquisition (Lawrence & Ochman, *PNAS* 1998).
- Contigs 5 and 6 (< 600 bp) contain no Prokka-called CDS and were excluded from downstream HMM scans.
