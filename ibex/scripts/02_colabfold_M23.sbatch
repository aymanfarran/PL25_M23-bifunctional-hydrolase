#!/bin/bash
#SBATCH --job-name=cf_M23
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|v100        # GPU constraint — adjust if needed
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/ibex/scratch/%u/non_endolysin/logs/cf_M23_%j.out
#SBATCH --error=/ibex/scratch/%u/non_endolysin/logs/cf_M23_%j.err

# ============================================================
# IBEX job: ColabFold structure prediction of PL25 M23 peptidase
# Predicts BOTH the full-length and the mature (no signal peptide) form
# ============================================================

set -euo pipefail

module load colabfold        # Adjust to actual module name on IBEX
# Or: source /ibex/sw/colabfold/setup.sh

WORK=/ibex/scratch/$USER/non_endolysin
IN=$WORK/inputs
OUT=$WORK/outputs/02_colabfold
mkdir -p "$OUT"

# Full-length M23
colabfold_batch \
    "$IN/PL25_M23.faa" \
    "$OUT/PL25_M23_full" \
    --num-models 5 \
    --num-recycle 3 \
    --model-type alphafold2_ptm \
    --amber \
    --templates

# Mature M23 (signal peptide removed) — more representative of secreted form
colabfold_batch \
    "$IN/PL25_M23_mature.faa" \
    "$OUT/PL25_M23_mature" \
    --num-models 5 \
    --num-recycle 3 \
    --model-type alphafold2_ptm \
    --amber \
    --templates

echo "Done. Top-ranked structures:"
ls -la "$OUT/PL25_M23_full"/*_rank_001*.pdb 2>/dev/null
ls -la "$OUT/PL25_M23_mature"/*_rank_001*.pdb 2>/dev/null
