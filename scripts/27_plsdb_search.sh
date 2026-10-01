#!/bin/bash
# Compare pPL25-1 (contig_11) against the PLSDB plasmid panel by two
# complementary methods:
#   (1) Mash sketch distance — fast genome-wide k-mer screen
#   (2) skani ANI            — sensitive ANI + alignment fraction
#
# Reference release used in the manuscript:
#   PLSDB 2025 (Schmartz et al.) — 72,556 closed plasmid sequences,
#   Mash sketch + per-plasmid FASTA, accessed 2026-06-18, from
#   https://ccb-microbe.cs.uni-saarland.de/plsdb/
#
# The database itself is NOT deposited in this repository: it is several GB and
# is redistributed by its authors. Download the release above and unpack it to
# $META (default: results/24_plasmid_db_search/plsdb_meta/) before running.
#
# Usage:  bash scripts/27_plsdb_search.sh
# Env:    MASH, SKANI, PYTHON override the executables; META overrides the
#         unpacked database location.
# Needs:  mash v2.3, skani v0.3.2, Biopython (only if the release ships one
#         concatenated FASTA rather than a fastas/ directory).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$ROOT/results/24_plasmid_db_search"
META="${META:-$WORK/plsdb_meta}"
QUERY="$WORK/pPL25-1.fna"

MASH="${MASH:-mash}"
SKANI="${SKANI:-skani}"
PYTHON="${PYTHON:-python3}"

[[ -f "$QUERY" ]] || { echo "query not found: $QUERY" >&2; exit 1; }
[[ -d "$META"  ]] || { echo "PLSDB release not found at $META — see header" >&2; exit 1; }

echo "[step 1] Mash sketch distance - pPL25-1 vs PLSDB"

MSH=$(find "$META" -maxdepth 3 -name "*.msh" -size +50M 2>/dev/null | head -1)
if [[ -z "$MSH" ]]; then
  echo "  WARN: PLSDB Mash sketch not found in $META"
else
  echo "  using sketch: $MSH"
  "$MASH" dist -p 4 "$MSH" "$QUERY" \
    | sort -k3,3g \
    > "$WORK/mash_dist_all.tsv"
  head -20 "$WORK/mash_dist_all.tsv" | tee "$WORK/mash_top20.tsv"
fi

echo ""
echo "[step 2] skani ANI - pPL25-1 vs PLSDB"

# Some PLSDB releases ship a fastas/ directory, others one concatenated FASTA.
PLSDB_FASTA_DIR=$(find "$META" -maxdepth 3 -type d -name "fastas" 2>/dev/null | head -1)
if [[ -z "$PLSDB_FASTA_DIR" ]]; then
  PLSDB_FASTA_FILE=$(find "$META" -maxdepth 3 -name "*.fasta" -size +500M 2>/dev/null | head -1)
  if [[ -z "$PLSDB_FASTA_FILE" ]]; then
    echo "  WARN: PLSDB FASTA(s) not found in $META; skani step skipped"
    exit 0
  fi
  echo "  splitting concatenated FASTA..."
  PLSDB_FASTA_DIR="$WORK/plsdb_fasta_split"
  mkdir -p "$PLSDB_FASTA_DIR"
  PLSDB_FASTA_FILE="$PLSDB_FASTA_FILE" PLSDB_FASTA_DIR="$PLSDB_FASTA_DIR" \
  "$PYTHON" - <<'PY'
from Bio import SeqIO
import os
src = os.environ["PLSDB_FASTA_FILE"]
dst = os.environ["PLSDB_FASTA_DIR"]
n = 0
for r in SeqIO.parse(src, "fasta"):
    SeqIO.write([r], os.path.join(dst, f"{r.id}.fasta"), "fasta")
    n += 1
print(f"  wrote {n} fasta files")
PY
fi

LST="$WORK/plsdb_fasta_list.txt"
find "$PLSDB_FASTA_DIR" -name "*.fasta" -o -name "*.fna" > "$LST"
echo "  $(wc -l < "$LST") plasmid fastas"

"$SKANI" search "$QUERY" --rl "$LST" \
  --min-af 30 -o "$WORK/skani_results.tsv" -t 4 2>&1 | tail -10
head -20 "$WORK/skani_results.tsv"

echo ""
echo "[done] outputs in $WORK"
