# IBEX batch scripts

SLURM scripts for the compute-heavy steps of the PL25_M23 analysis: the
remote BLASTp searches, the AlphaFold2 model, the geNomad scan and the
initial M23 phylogeny. Everything downstream of these runs locally from
`../scripts/`.

## Scripts

| Script | What it does | Inputs | Outputs | Walltime | GPU |
|---|---|---|---|---|---|
| `01_blast_nr.sbatch` | BLASTp of three PL25 queries against NCBI nr | `PL25_M23.faa`, `PL25_00076_unknown.faa`, `contig_11_proteome.faa` | `outputs/01_blast_nr/*_vs_nr.tsv` | ~3 h | no |
| `02_colabfold_M23.sbatch` | AlphaFold2 structure of PL25_M23, full length and residues 38–370 | `PL25_M23.faa`, `PL25_M23_mature.faa` | `outputs/02_colabfold/PL25_M23_{full,mature}/` | ~3 h | **yes** |
| `03_genomad_PL25.sbatch` | geNomad plasmid/virus classification of the whole assembly | `PL25.fasta` | `outputs/03_genomad/` | ~1 h | no |
| `04_iqtree_M23_phylogeny.sbatch` | MAFFT → trimAl → IQ-TREE, 1000 UFBoot | `M23_seedset.faa` (built from script 01) | `outputs/04_M23_phylogeny/M23_iqtree.*` | ~1 h | no |
| `05_iceberg_blast.sbatch` | BLASTn of contig_11 against ICEberg2 — **exploratory, not used for the conclusions in the manuscript** | `contig_11.fna` | `outputs/05_iceberg/` | ~30 m | no |

The manuscript models the mature form (residues 38–370) produced by script 02;
the full-length model is retained for reference only.

## Software

| Tool | Module |
|---|---|
| BLAST+ | `blast/2.16.0` |
| MAFFT | `mafft/7.526` |
| trimAl | `trimal/1.4.1` |
| IQ-TREE | `iqtree/2.3.4` |
| ColabFold | `colabfold` |
| conda | `miniforge` (for geNomad) |

> **Note on the phylogeny.** Script 04 produces the initial M23 tree under
> IQ-TREE 2.3.4 with automatic model selection (`-m MFP`). The tree reported
> in the manuscript was rebuilt locally under IQ-TREE v3.1.1 with the
> Q.pfam+G4 model; see the phylogeny section of the top-level `README.md`.

## Order

Script 01 must finish first: its BLASTp output supplies the seed set for the
phylogeny. Scripts 02, 03 and 05 are independent and can run in parallel with
it.

## Setup

On IBEX:

```bash
WORK=/ibex/scratch/$USER/non_endolysin
mkdir -p $WORK/{inputs,outputs,logs,scripts}
```

From a local terminal:

```bash
scp -r ibex/inputs  $USER@ilogin:/ibex/scratch/$USER/non_endolysin/
scp -r ibex/scripts $USER@ilogin:/ibex/scratch/$USER/non_endolysin/
```

Before the first submission, confirm the module names resolve and set the
database paths inside the scripts — `NR_DB` in script 01, the geNomad
database in script 03, and the ICEberg2 download in script 05:

```bash
module avail blast mafft iqtree colabfold miniforge
genomad download-database /ibex/scratch/$USER/genomad_db
```

## Running

```bash
cd $WORK/scripts

sbatch 01_blast_nr.sbatch
sbatch 02_colabfold_M23.sbatch
sbatch 03_genomad_PL25.sbatch
sbatch 05_iceberg_blast.sbatch

squeue -u $USER
```

Once script 01 has finished, build the phylogeny seed set from its output:

- the top 80 high-scoring unique PL25_M23 homologues from `PL25_M23_vs_nr.tsv`
- the five characterised reference M23 enzymes — lysostaphin `P10547`,
  LytM `O33599`, ALE-1 `O05156`, EnpA `E1V3I0`, zoocin A `O54309`
- the CwlP catalytic core (residues 1373–1686 of `O31976`)
- 14 phage-associated M23 lysins

giving 101 sequences. Save as
`$WORK/outputs/04_M23_phylogeny/M23_seedset.faa`, then:

```bash
sbatch 04_iqtree_M23_phylogeny.sbatch
```

CwlT (`P96645`) is deliberately excluded: its C-terminal peptidase belongs to
the NlpC/P60 family (PF00877), not Peptidase_M23 (PF01551). It is used only as
a structural comparator, in `../scripts/05cc_M23_3way_alignment.py`.

## Retrieving results

```bash
rsync -avz $USER@ilogin:/ibex/scratch/$USER/non_endolysin/outputs/ ibex/outputs/
```
