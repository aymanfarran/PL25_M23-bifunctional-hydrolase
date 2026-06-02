#!/bin/bash
#SBATCH --job-name=blast_nr_M23
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=/ibex/scratch/%u/non_endolysin/logs/blast_nr_%j.out
#SBATCH --error=/ibex/scratch/%u/non_endolysin/logs/blast_nr_%j.err

# ============================================================
# IBEX job: BLASTp 3 PL25 queries vs NCBI nr database
# Targets:
#   - PL25_M23 (PL25_00078)       — find closest M23 homologues for phylogeny
#   - PL25_00076_unknown          — identify Rep vs Relaxase
#   - contig_11.fna  (nucleotide BLAST-X vs nr to confirm gene context)
# ============================================================

set -euo pipefail

module load blast/2.16.0   # Adjust to whatever's on IBEX

# Path to local nr (check actual location on IBEX)
NR_DB="/ibex/reference/KSL/blast_db/nr"          # ← VERIFY THIS PATH on IBEX
# Alternate: NR_DB=/scratch/ksl/refseq/nr

WORK=/ibex/scratch/$USER/non_endolysin
IN=$WORK/inputs
OUT=$WORK/outputs/01_blast_nr
mkdir -p "$OUT" "$WORK/logs"

# ── 1. M23 protein vs nr (for phylogeny seed set) ────────────────────
blastp -query "$IN/PL25_M23.faa" \
       -db "$NR_DB" \
       -out "$OUT/PL25_M23_vs_nr.tsv" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids" \
       -evalue 1e-10 \
       -max_target_seqs 500 \
       -num_threads 16

# Also get full alignment for top hits
blastp -query "$IN/PL25_M23.faa" \
       -db "$NR_DB" \
       -out "$OUT/PL25_M23_vs_nr.txt" \
       -outfmt 0 \
       -evalue 1e-10 \
       -max_target_seqs 100 \
       -num_threads 16

# ── 2. PL25_00076 (unknown Rep/Relaxase) vs nr ────────────────────────
blastp -query "$IN/PL25_00076_unknown.faa" \
       -db "$NR_DB" \
       -out "$OUT/PL25_00076_vs_nr.tsv" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids" \
       -evalue 1e-5 \
       -max_target_seqs 200 \
       -num_threads 16
blastp -query "$IN/PL25_00076_unknown.faa" \
       -db "$NR_DB" \
       -out "$OUT/PL25_00076_vs_nr.txt" \
       -outfmt 0 \
       -evalue 1e-5 \
       -max_target_seqs 50 \
       -num_threads 16

# ── 3. contig_11 nucleotide vs nr (BLASTX, find related ICEs/plasmids) ───
blastx -query "$IN/contig_11.fna" \
       -db "$NR_DB" \
       -out "$OUT/contig_11_blastx_nr.tsv" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids" \
       -evalue 1e-10 \
       -max_target_seqs 50 \
       -num_threads 16 \
       -query_gencode 11

echo "Done. Results in $OUT/"
ls -la "$OUT/"
