#!/bin/bash
#SBATCH --job-name=iceberg_blast
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/ibex/scratch/%u/non_endolysin/logs/iceberg_%j.out
#SBATCH --error=/ibex/scratch/%u/non_endolysin/logs/iceberg_%j.err

# ============================================================
# IBEX job: BLAST contig_11 against ICEberg2 database
# Goal: find closest characterised ICE/plasmid for context
# ============================================================

set -euo pipefail

module load blast/2.16.0

WORK=/ibex/scratch/$USER/non_endolysin
IN=$WORK/inputs
DB=/ibex/scratch/$USER/iceberg2_db    # Pre-download from https://bioinfo-mml.sjtu.edu.cn/ICEberg2/

OUT=$WORK/outputs/05_iceberg
mkdir -p "$OUT"

# Pre-build BLAST db from ICEberg FASTA if not done
if [ ! -f "$DB.nhr" ]; then
    makeblastdb -in "$DB.fna" -dbtype nucl -out "$DB"
fi

# contig_11 nucleotide vs ICEberg
blastn -query "$IN/contig_11.fna" \
       -db "$DB" \
       -out "$OUT/contig_11_vs_iceberg.tsv" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
       -evalue 1e-10 \
       -max_target_seqs 50 \
       -num_threads 8

# contig_11 proteome vs ICEberg proteins
blastp -query "$IN/contig_11_proteome.faa" \
       -db "$DB.prot" \
       -out "$OUT/contig_11_proteome_vs_iceberg.tsv" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
       -evalue 1e-5 \
       -max_target_seqs 5 \
       -num_threads 8

echo "Done. Top hits:"
sort -k12,12 -gr "$OUT/contig_11_vs_iceberg.tsv" | head -10
