# Supplementary Methods S7 — Closed-genome PL25_M23 homologue survey and focused phylogeny

**Closed-genome survey.** For each of the 250 BLASTp hits, source assemblies were retrieved from the NCBI Identical Protein Group (IPG) report (`efetch -db protein -rettype ipg -retmode xml`). Assemblies were retained only if their `assembly_level` was `Complete Genome` or `Chromosome` (`esummary -db assembly`) and their source genus matched a curated list of halotolerant Bacillaceae: *Halobacillus*, *Lentibacillus*, *Virgibacillus*, *Gracilibacillus*, *Oceanobacillus*, *Cytobacillus*, *Ornithinibacillus*, *Salinicoccus*, *Salimicrobium*, *Aquibacillus*, *Pontibacillus*, *Halolactibacillus*, *Salibacillus*, *Salinibacillus*, *Salirhabdus*. Replicon type was assigned by searching the source nucleotide FASTA defline for *"plasmid"* or *"chromosome"*. Fourteen assemblies passed both filters (Supplementary Table S5).

**Focused phylogeny.** Full-length protein sequences of the 14 homologues, PL25_M23, and lysostaphin (UniProt P10547) were aligned with MAFFT v7.526 (`--auto`; 1,071 columns), trimmed with trimAl v1.5.rev1 (`-gappyout`; 332 informative columns retained), and used to infer a maximum-likelihood phylogeny with IQ-TREE v3.1.1 under WAG+G4 selected by ModelFinder under BIC (Whelan & Goldman, 2001). Branch support was estimated from 1,000 ultrafast bootstrap (UFBoot; Hoang et al., 2018) and 1,000 SH-aLRT replicates (Anisimova et al., 2011). The tree was rooted on lysostaphin and rendered with Biopython 1.87 Bio.Phylo. Scripts and intermediate files are provided at `results/22_closed_genome_search/` in the project Zenodo archive.

## References (new in this section)

MAFFT, trimAl, IQ-TREE, ModelFinder, Biopython and NCBI database references are listed in the main bibliography. The three new references introduced by this analysis are:

- Anisimova M, Gil M, Dufayard JF, Dessimoz C, Gascuel O. Survey of branch support methods demonstrates accuracy, sensitivity, and reproducibility of fast likelihood-based approximation schemes. *Syst Biol*. 2011;60:685–699.
- Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. UFBoot2: Improving the ultrafast bootstrap approximation. *Mol Biol Evol*. 2018;35:518–522.
- Whelan S, Goldman N. A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. *Mol Biol Evol*. 2001;18:691–699.
