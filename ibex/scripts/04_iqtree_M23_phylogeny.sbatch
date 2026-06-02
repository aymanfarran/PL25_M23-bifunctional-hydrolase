#!/bin/bash
#SBATCH --job-name=iqtree_M23
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=/ibex/scratch/%u/non_endolysin/logs/iqtree_M23_%j.out
#SBATCH --error=/ibex/scratch/%u/non_endolysin/logs/iqtree_M23_%j.err

# ============================================================
# IBEX job: M23 phylogeny — PL25_M23 + UniRef90 homologues
#   1. MAFFT alignment
#   2. trimAl  cleanup
#   3. IQ-TREE2 with ModelFinder + 1,000 UFBoot + 1,000 SH-aLRT
# ============================================================

set -euo pipefail

module load mafft/7.526
module load trimal/1.4.1
module load iqtree/2.3.4

WORK=/ibex/scratch/$USER/non_endolysin
OUT=$WORK/outputs/04_M23_phylogeny
mkdir -p "$OUT"

cd "$OUT"

# Input: combined FASTA of PL25_M23 + closest homologues
# This file is produced AFTER BLAST step 01 — combine top 100 hits + reference seqs
INPUT_FASTA="$OUT/M23_seedset.faa"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: $INPUT_FASTA not found."
    echo "Build it from blastp results (step 01) + manually-added references:"
    echo "  - Lysostaphin (P10547)"
    echo "  - LytM (Q8NXI8)"
    echo "  - ALE-1 (Q47728)"
    echo "  - CwlT (O32008)"
    echo "  - EnpA (E1V3I0)"
    echo "  - Zoocin A (O54309)"
    exit 1
fi

# 1. Align
mafft --auto --thread 16 "$INPUT_FASTA" > M23_aln.fasta

# 2. Trim
trimal -in M23_aln.fasta -out M23_aln_trim.fasta -gappyout

# 3. ML tree with ModelFinder + bootstrap
iqtree2 \
    -s M23_aln_trim.fasta \
    -m MFP \
    -bb 1000 -alrt 1000 \
    -nt 16 \
    -pre M23_iqtree

echo "Done. Tree: $OUT/M23_iqtree.treefile"
