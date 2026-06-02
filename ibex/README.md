# IBEX scripts — non_endolysin paper (M23 only)

## Overview
Five SLURM batch scripts for analyses too heavy to run locally. Run **after** you log in to IBEX.

| Script | What it does | Walltime | GPU? |
|---|---|---|---|
| `01_blast_nr.sbatch`         | BLASTp 3 queries vs NCBI nr  | ~3 h  | no  |
| `02_colabfold_M23.sbatch`    | AlphaFold2 structure of M23  | ~3 h  | **yes** |
| `03_genomad_PL25.sbatch`     | geNomad full-genome scan     | ~1 h  | no  |
| `04_iqtree_M23_phylogeny.sbatch` | IQ-TREE2 1000 UFBoot tree | ~1 h  | no  |
| `05_iceberg_blast.sbatch`    | BLAST contig_11 vs ICEberg2  | ~30 m | no  |

## Step 0 — Set up the IBEX scratch tree (one-time)

After SSH'ing into IBEX:

```bash
# Replace $USER with your IBEX username
WORK=/ibex/scratch/$USER/non_endolysin
mkdir -p $WORK/{inputs,outputs,logs,scripts}

# Copy inputs and scripts from your Mac:
# Run THIS from your Mac terminal (not IBEX):
#   scp -r ~/prophage-endolysin-pipeline/non_endolysin_paper/ibex/inputs  $USER@ilogin:/ibex/scratch/$USER/non_endolysin/
#   scp -r ~/prophage-endolysin-pipeline/non_endolysin_paper/ibex/scripts $USER@ilogin:/ibex/scratch/$USER/non_endolysin/

cd $WORK/scripts
```

## Step 1 — Run BLAST vs nr FIRST (depends on this)

```bash
# Verify NR_DB path inside the script — edit if needed
sbatch 01_blast_nr.sbatch
squeue -u $USER
```

Output: `$WORK/outputs/01_blast_nr/`
Critical files:
- `PL25_M23_vs_nr.tsv`  — feeds phylogeny (script 04)
- `PL25_00076_vs_nr.tsv` — resolves Rep vs Relaxase question

## Step 2 — Submit the other jobs in parallel

Once script 01 completes (or in parallel for the ones that don't depend on it):

```bash
# Independent — can submit in parallel
sbatch 02_colabfold_M23.sbatch
sbatch 03_genomad_PL25.sbatch
sbatch 05_iceberg_blast.sbatch

# DEPENDS ON script 01: phylogeny needs the BLAST seed set first
# After script 01 finishes, manually build M23_seedset.faa from top BLAST hits:
#   - Take top 100 unique hits from PL25_M23_vs_nr.tsv
#   - Pull their UniRef90 sequences (use efetch or download)
#   - Add reference enzymes (lysostaphin P10547, LytM Q8NXI8, ALE-1 Q47728, CwlT O32008, EnpA E1V3I0)
#   - Save to $WORK/outputs/04_M23_phylogeny/M23_seedset.faa
sbatch 04_iqtree_M23_phylogeny.sbatch
```

## Pre-flight checks BEFORE submitting

Each script has paths/modules to verify on IBEX:

1. **`01_blast_nr.sbatch`** — check `NR_DB` path is correct
2. **`02_colabfold_M23.sbatch`** — verify `module load colabfold` syntax; check GPU partition name
3. **`03_genomad_PL25.sbatch`** — pre-download geNomad db: `genomad download-database /ibex/scratch/$USER/genomad_db`
4. **`05_iceberg_blast.sbatch`** — pre-download ICEberg2 from https://bioinfo-mml.sjtu.edu.cn/ICEberg2/

Run this check once on IBEX:
```bash
module avail blast
module avail mafft
module avail iqtree
module avail colabfold     # or: ls /ibex/sw/colabfold
module avail miniforge     # or whichever conda installer
```

## Quick smoke tests (1-minute jobs)

To verify environment before running the real jobs:

```bash
# Smoke test BLAST
sbatch --time=10 --wrap="module load blast/2.16.0 && blastp -version && which blastp"

# Smoke test conda for geNomad
sbatch --time=10 --wrap="conda activate genomad && genomad --version"

# Smoke test IQ-TREE
sbatch --time=10 --wrap="module load iqtree/2.3.4 && iqtree2 -version"
```

## When jobs finish — pull results back to Mac

```bash
# Run on your Mac:
LOCAL=~/prophage-endolysin-pipeline/non_endolysin_paper/ibex/outputs
mkdir -p $LOCAL
rsync -avz $USER@ilogin:/ibex/scratch/$USER/non_endolysin/outputs/ $LOCAL/
```

After this, downstream analyses (annotation, HMM scans, geNomad, IQ-TREE, ColabFold) can be re-run locally or on IBEX from `non_endolysin_paper/scripts/`.

## Note on M15
M15 was dropped from this paper (not actually a PL25 protein). Only M23 is being characterised.
