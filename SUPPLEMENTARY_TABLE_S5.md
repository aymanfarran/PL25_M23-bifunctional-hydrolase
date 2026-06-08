# Supplementary Table S5

**PL25_M23 homologues in closed NCBI Bacillaceae assemblies: replicon context, protein hit, and BLASTp statistics.**

Closed (Complete Genome or Chromosome-level) NCBI assemblies of halotolerant Bacillaceae carrying a PL25_M23 BLAST homologue were identified by querying the NCBI Identical Protein Group (IPG) report for each of the 250 BLASTp hits and retaining only assemblies whose source genus is recognised as halotolerant Bacillaceae (see Supplementary Methods S7). Replicon type was assigned from the title of each source nucleotide record. Plasmid-encoded entries are marked **★**; chromosomal entries carry no marker. Eleven of fourteen (79 %) homologues are plasmid-encoded across five of the seven represented genera.

| # | Organism (strain) | Assembly accession | Replicon | Plasmid name | Source nucleotide | Protein hit | Length (aa) | % identity to PL25_M23 | BLASTp E-value |
|---|---|---|:---:|---|---|---|---:|---:|---:|
| 1 | *Cytobacillus oceanisediminis* 2691 | GCF_000294775.2 | ★ plasmid | pBO1 | NZ_CP015507.1 | WP_019380455.1 | 703 | 44.8 % | 4.7 × 10⁻⁶⁵ |
| 2 | *Halobacillus litoralis* ERB 031 | GCF_004101865.1 | ★ plasmid | pLDW-31 | NZ_CP026119.1 | WP_241656025.1 | 620 | 43.8 % | 1.8 × 10⁻⁶¹ |
| 3 | *Virgibacillus halodenitrificans* Bac324 | GCF_003667805.1 | ★ plasmid | unnamed | NZ_CP033050.1 | WP_121616750.1 | 758 | 48.4 % | 3.8 × 10⁻⁶⁵ |
| 4 | *Gracilibacillus boraciitolerans* JCM 21714 | GCF_055385745.1 | ★ plasmid | pJCM21714 | NZ_AP043914.1 | WP_464380148.1 | 363 | 50.0 % | 1.1 × 10⁻¹⁰⁶ |
| 5 | *Salinicoccus jeotgali* JCM 16981 | GCF_055379195.1 | ★ plasmid | pJCM16981-2 | NZ_AP043832.1 | WP_344700570.1 | 426 | 35.0 % | 6.7 × 10⁻⁶¹ |
| 6 | *Cytobacillus firmus* M7 | GCF_025732135.1 | ★ plasmid | p1 | NZ_CP107028.1 | WP_263600221.1 | 707 | 44.5 % | 2.7 × 10⁻⁶⁴ |
| 7 | *Cytobacillus firmus* CK19 | GCF_023823815.1 | ★ plasmid | unnamed1 | NZ_CP085391.1 | WP_252222728.1 | 707 | 44.5 % | 2.9 × 10⁻⁶⁴ |
| 8 | *Cytobacillus firmus* SGAir0285 | GCF_043540845.1 | ★ plasmid | pSGAir0285_3 | NZ_CP027985.1 | WP_009335593.1 | 627 | 43.4 % | 1.7 × 10⁻⁶¹ |
| 9 | *Cytobacillus oceanisediminis* CK22 | GCF_023823835.1 | ★ plasmid | unnamed1 | NZ_CP085393.1 | WP_252246573.1 | 707 | 44.5 % | 3.6 × 10⁻⁶⁴ |
| 10 | *Cytobacillus oceanisediminis* YPW-V2 | GCF_014883935.1 | ★ plasmid | unnamed1 | NZ_CP062791.1 | WP_192909257.1 | 627 | 43.8 % | 2.4 × 10⁻⁶¹ |
| 11 | *Cytobacillus pseudoceanisediminis* | GCF_023516215.1 | ★ plasmid | unnamed | NZ_CP097350.1 | WP_251264864.1 | 707 | 44.5 % | 2.3 × 10⁻⁶³ |
| 12 | *Lentibacillus* sp. CBA3610 | GCF_013373365.1 | chromosome | — | NZ_CP035925.1 | WP_176447356.1 | 361 | 52.9 % | 8.2 × 10⁻¹²⁸ |
| 13 | *Lentibacillus amyloliquefaciens* LAM0015 | GCF_001307805.1 | chromosome | — | NZ_CP013862.1 | WP_068441797.1 | 358 | 51.1 % | 6.7 × 10⁻¹²³ |
| 14 | *Virgibacillus ainsalahensis* JCM 30907 | GCF_055380755.1 | chromosome | — | NZ_AP043943.1 | WP_461177648.1 | 360 | 53.0 % | 6.9 × 10⁻¹²⁶ |

## Notes

1. **Replicon assignment** was made by parsing the title of each source nucleotide record (NCBI esummary `db=nuccore`); entries with "plasmid" in the title were classified as plasmid-encoded, the remainder as chromosomal (titles ending "chromosome, complete genome").
2. **Multi-domain plasmid-encoded homologues**: plasmid-encoded homologues range from 363 aa to 758 aa, with most (n = 9) in the 620–758 aa range, compared with chromosomal homologues (358–361 aa) and PL25_M23 itself (370 aa). The larger plasmid-encoded variants are likely multi-domain fusion proteins; the identity of their additional N-/C-terminal domains was not investigated in this study.
3. **% identity and E-value** are taken from the original BLASTp search of PL25_M23 (38–370 aa) against NCBI nr (`-evalue 1e-5 -max_target_seqs 250`); see main Methods. Identity is to the full-length BLAST alignment region, not to the PL25_M23 M23 catalytic domain alone.
4. **Halotolerant Bacillaceae genera** considered: *Halobacillus*, *Lentibacillus*, *Virgibacillus*, *Gracilibacillus*, *Oceanobacillus*, *Cytobacillus*, *Ornithinibacillus*, *Salinicoccus*, *Salimicrobium*, *Aquibacillus*, *Pontibacillus*, *Halolactibacillus*, *Salibacillus*, *Salinibacillus*, *Salirhabdus*. Only the seven genera represented in the closed-assembly hit list are shown above.
