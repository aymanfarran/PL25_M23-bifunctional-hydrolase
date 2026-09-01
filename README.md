# PL25_M23 — a plasmid-candidate-encoded bifunctional SLT/M23 hydrolase from *Virgibacillus salarius* PL25 with activity against Gram-positive bacteria

This repository contains the code, source data, intermediate files, structural models and figure-generation scripts supporting:

> [Author list]. *A divergent plasmid-encoded SLT/M23 cell-wall hydrolase from halophilic Virgibacillus salarius PL25 with antibacterial activity against Gram-positive pathogens.* [Journal], [Year]. DOI: [paper DOI when assigned].

---

## Quick reference

| Item | Identifier |
|---|---|
| BioProject | [PRJNA1473924](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1473924) |
| BioSample | [SAMN60541943](https://www.ncbi.nlm.nih.gov/biosample/SAMN60541943) |
| WGS assembly | [JBYTNQ000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBYTNQ000000000) (version `JBYTNQ010000000`) |
| Raw Nanopore reads (SRA) | [SRR35929825](https://www.ncbi.nlm.nih.gov/sra/SRR35929825) |
| Locus tag prefix | `AC3VG8` |
| Strain | *Virgibacillus salarius* PL25 (TaxID 447199) |
| Source | Hypersaline-lake sediment, Al Lith, Saudi Arabia (2023) |
| Sequencing | Oxford Nanopore PromethION; basecalled with Guppy (high-accuracy) |
| Assembly | Flye v2.9.3-b1797, `--nano-raw` |

---

## Repository structure

```
.
├── LICENSE                          MIT license
├── README.md                        This file
├── .gitignore
│
├── scripts/                         All analysis and plotting scripts
│   ├── 04_realign_ALE1.pml                 PyMOL — ALE-1 M23 superposition
│   ├── 04c_render_M23_panels.py            PyMOL — Figure 3b structure views
│   ├── 04c_combine_M23_panel.py            4-view 3D structure assembly
│   ├── 05_active_site_panel.py             PyMOL — Zn²⁺-site close-ups
│   ├── 05b_combine_panels.py               Active-site panel assembly
│   ├── 05cc_M23_3way_alignment.py          PyMOL — Figure 4 panel rendering
│   ├── 05cc_combine.py                     Figure 4 — three-way alignment assembly
│   ├── 06c_M23_bifunctional_figure_v3.py   Figure 3a — pLDDT + domain architecture
│   ├── 08_M23_tree_v4.py                   Supp. Figure 1 — full 101-taxon phylogeny
│   ├── 08_render_surfaces.py               PyMOL — Figure 5 surface rendering
│   ├── 08_combine_electrostatic.py         Figure 5 — electrostatic comparison assembly
│   ├── 09_circular_plasmid_map.py          Figure 1a — circular plasmid map
│   ├── 10_antibacterial_activity.py        Primary antibacterial assay stats
│   ├── 11_orit_motif_scan.py               oriT consensus-motif screen (Supp. Methods S4)
│   ├── 22_closed_genome_search.py          Closed-genome M23 homologue census
│   ├── 22_render_F3_plasmid_lineage.py     Supp. Figure 2 — plasmid-clade phylogeny
│   ├── 39_compact_M23_targeted.py          Figure 2 — compact targeted M23 phylogeny
│   ├── 40_activity_export_raw.py           Bench workbook → tidy raw-data TSVs
│   └── 41_activity_panels.py               Figures 6 and 7 — activity + biochemistry
│
│   (figures/ is not tracked in this repository — see "Figures" below)
│
├── results/
│   ├── 01_annotation/prokka/        Prokka genome annotation (GBK, GFF, FAA, FFN)
│   ├── 03_M23_ICE_locus/            contig_11 sequence and CDS
│   ├── 04_HMM_scans/                Custom 12-profile Pfam HMM panel + scan output
│   │   └── pfam_db/                  PF00263, PF00589, PF01551, PF02534, PF03135,
│   │                                 PF03389, PF03432, PF03524, PF13245,
│   │                                 PF13407, PF13476, PF13701
│   ├── 05_blast_nr/                 Remote BLASTp results for PL25_00076–00148
│   ├── 06_genomad/                  geNomad v1.7.0 classification
│   ├── 07_M23_phylogeny/            MAFFT alignment, trimAl, IQ-TREE v4 tree (101 seqs)
│   ├── 08_colabfold/                AlphaFold2 model of PL25_M23 (residues 38–370)
│   ├── 09_gc_codon/                 GC / GC1 / GC2 / GC3 per contig
│   ├── 10_ncbi_plasmids/            NCBI WGS GenBank records for synteny
│   ├── 11_interproscan/             InterProScan v5.61-93.0 domain output
│   ├── 12_signal_peptide/           SignalP 6.0h output
│   ├── 13_M23_refs/                 Reference structures (PDB + AF2 models)
│   ├── 14_cwlp_structure/           AlphaFold2 model of CwlP catalytic core
│   ├── 14_synteny/                  clinker synteny analysis
│   ├── 15_mob_suite/                MOB-suite v3.1.9 typing
│   ├── 16_electrostatic/            pI / charge / charge-density JSON
│   ├── 17_conjscan/                 MacSyFinder + CONJScan v2.1.0 output
│   ├── 19_circularity/              FASTA-level circularity check (nucmer + oriT scan output)
│   └── 20_activity/                 Antibacterial assay raw data + stats
│
└── ibex/
    ├── README.md                    IBEX usage notes
    ├── scripts/                     SLURM job scripts
    └── inputs/                      Input FASTAs sent to IBEX jobs
```

**Note:** the manuscript prose, the rendered figures and the supplementary tables are all maintained outside this repository and are distributed with the paper. The deposit holds the **code, source data, intermediate files, structural models and figure-generation scripts** needed to reproduce the analyses and rebuild the figures.

---

## Figures

Rendered figures are **not tracked in this repository**. They are distributed with
the manuscript and its supplementary material, and every one of them can be
regenerated from the scripts and result files kept here.

Running any script below writes its output into a local `figures/` directory,
which is git-ignored.

## How to reproduce a single figure

Most figures are one script, or a PyMOL render followed by a matplotlib
assembly step. Figures 1 and 3 are composited by hand from the panels their
scripts produce; Figure 6 panel b is a photograph of the drop plates and has
no script.

```bash
# Figure 1a — circular plasmid map
python scripts/09_circular_plasmid_map.py

# Figure 2 — compact targeted M23 phylogeny
python scripts/39_compact_M23_targeted.py

# Figure 3a — pLDDT + domain architecture
python scripts/06c_M23_bifunctional_figure_v3.py

# Figure 3b — structure views (requires PyMOL Open-Source)
pymol -cq scripts/04c_render_M23_panels.py

# Figure 4 — three-way structural alignment (requires PyMOL Open-Source)
pymol -cq scripts/05cc_M23_3way_alignment.py
python scripts/05cc_combine.py

# Figure 5 — electrostatic surface comparison (requires PyMOL Open-Source)
pymol -cq scripts/08_render_surfaces.py
python scripts/08_combine_electrostatic.py

# Figures 6 and 7 — activity assays and biochemical characterisation
python scripts/41_activity_panels.py
```

PyMOL steps must be run under the **open-source** build. The Incentive/trial
build stamps an evaluation watermark into the ray-traced PNGs, which the
assembly step would then bake into the composite figure.

The activity figures read the raw-data TSVs in `results/20_activity/`, which
are tracked here. To regenerate those TSVs from the bench workbook:

```bash
python scripts/40_activity_export_raw.py path/to/data.xlsx
```

The CwlT-free M23 phylogeny can be rebuilt from scratch with:

```bash
cd results/07_M23_phylogeny
mafft --auto M23_seedset_v4_noCwlT.faa > M23_aln_v4.fasta
trimal -in M23_aln_v4.fasta -out M23_aln_v4_trim.fasta -gappyout
iqtree3 -s M23_aln_v4_trim.fasta -m Q.pfam+G4 -B 1000 -T AUTO --prefix M23_iqtree_v4
python ../../scripts/08_M23_tree_v4.py          # Supp. Figure 1, all 101 taxa
python ../../scripts/39_compact_M23_targeted.py # Figure 2, collapsed to 20 tips
```

---

## Software versions

The exact tool versions and parameters used in the analyses are listed in the manuscript Methods. Briefly:

| Tool | Version |
|---|---|
| Flye | 2.9.3-b1797 |
| Prokka | 1.15.6 (Prodigal 2.6.3) |
| geNomad | 1.7.0 |
| MOB-suite | 3.1.9 |
| MacSyFinder + CONJScan | 2.1.6 + 2.1.0 |
| HMMER | 3.4 |
| BLAST+ | 2.16.0 |
| InterProScan | 5.61-93.0 |
| SignalP | 6.0h |
| MAFFT | 7.526 |
| trimAl | 1.5.rev1 |
| IQ-TREE | 3.1.1 |
| ColabFold (AlphaFold2 ptm) | 1.5.5 |
| PyMOL (Open-Source) | 3.1.6.1 |
| clinker | 0.0.32 |
| Biopython | 1.83 |
| Python | 3.11 |

---

## Data dependencies (databases consulted)

| Database | Version / access date |
|---|---|
| Pfam | release 37.0, downloaded 2026-05-14 |
| NCBI nr | BLASTp queried 2026-05-17 |
| UniProt Knowledgebase | accessions retrieved 2026-05-18 |
| RCSB PDB | 4QPB, 4ZYB, 6S6E, 5KVP downloaded 2026-05-17 |
| AlphaFold Protein Structure Database | AF-O05156-F1 v6 and AF-P96645-F1 v4 accessed 2026-05-22 |
| NCBI WGS (synteny contigs) | NZ_LDPV02000047.1, NZ_JARSFA010000044.1, NZ_JBCITH010000002.1 accessed 2026-05-17 |
| geNomad reference DB | v1.7 (used 2026-05-17) |

---

## License

- **Code (scripts, configs, notebooks):** MIT License (see `LICENSE`).
- **Data files, figures and supplementary tables:** CC-BY-4.0 — please cite the paper above when reusing.

---

## Citation

If you use any code, data, or figures from this repository, please cite the paper above and this archive (Zenodo DOI: *to be assigned upon first release*).
