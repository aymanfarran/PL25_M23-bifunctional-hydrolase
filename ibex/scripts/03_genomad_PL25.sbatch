#!/bin/bash
#SBATCH --job-name=genomad_PL25
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=/ibex/scratch/%u/non_endolysin/logs/genomad_%j.out
#SBATCH --error=/ibex/scratch/%u/non_endolysin/logs/genomad_%j.err

# ============================================================
# IBEX job: geNomad full-genome scan of PL25
# Goal: independent classification of contig_11 — plasmid vs ICE vs other
# geNomad replaces VirSorter2 + PlasmidFinder + ICEberg as a unified ML tool
# ============================================================

set -euo pipefail

module load miniforge        # Or however you load conda on IBEX
# conda activate genomad     # Pre-create env: conda create -n genomad -c conda-forge -c bioconda genomad

WORK=/ibex/scratch/$USER/non_endolysin
IN=$WORK/inputs
OUT=$WORK/outputs/03_genomad
DB=/ibex/scratch/$USER/genomad_db   # Pre-download: genomad download-database $DB

mkdir -p "$OUT"

# Run end-to-end with permissive thresholds to catch any divergent ICEs
genomad end-to-end \
    --cleanup \
    --splits 8 \
    --threads 16 \
    --min-virus-marker-enrichment 0 \
    --min-plasmid-marker-enrichment 0 \
    --min-virus-hallmarks 0 \
    --min-plasmid-hallmarks 0 \
    "$IN/PL25.fasta" \
    "$OUT" \
    "$DB"

echo "Done. Key outputs:"
ls "$OUT/PL25_summary"/ 2>/dev/null
