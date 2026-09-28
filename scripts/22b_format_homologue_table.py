#!/usr/bin/env python3
"""Format the closed-genome homologue census into the supplementary table.

Turns the machine-readable output of 22_closed_genome_search.py into the
presentation table used in the supplementary material: adds a stable rank and
the protein length, and renders percent identity and e-values for print. The
host label, including each strain designation, is taken from the census as
written, so the table and the phylogeny tips always name the same isolates.

Inputs (results/22_closed_genome_search/):
  • closed_M23_homologues.tsv   — census written by 22_closed_genome_search.py
  • closed_M23_homologues.faa   — full-length sequences, for length_aa

Output:
  • closed_M23_homologues_table.tsv

Row order is plasmid-encoded homologues first, then chromosomal, each block
sorted by descending percent identity to PL25_M23. This is a documented
ordering rather than the order NCBI happened to return the records in, so the
table is byte-identical on every rerun.

Usage:  python scripts/22b_format_homologue_table.py
Requires: Python >= 3.10, no third-party packages.
"""
from __future__ import annotations
import csv, sys
from pathlib import Path

ROOT   = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "results/22_closed_genome_search"
CENSUS = OUTDIR / "closed_M23_homologues.tsv"
FASTA  = OUTDIR / "closed_M23_homologues.faa"
TABLE  = OUTDIR / "closed_M23_homologues_table.tsv"

COLS = ["rank", "organism", "assembly", "replicon", "plasmid_name",
        "nucleotide_acc", "protein_acc", "length_aa", "pct_identity", "evalue"]

SUPERSCRIPT = str.maketrans("0123456789-", "⁰¹²³⁴"
                                           "⁵⁶⁷⁸⁹⁻")


def seq_lengths(path: Path) -> dict[str, int]:
    """Protein accession -> sequence length, from the census FASTA."""
    lengths, acc = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            acc = line[1:].split("|")[-1]
            lengths[acc] = 0
        elif acc:
            lengths[acc] += len(line.strip())
    return lengths


def fmt_evalue(ev: str) -> str:
    """'4.70e-65' -> '4.7 x 10^-65' with typographic superscripts."""
    try:
        mant, exp = f"{float(ev):.1e}".split("e")
    except ValueError:
        return ev
    return f"{mant} × 10{str(int(exp)).translate(SUPERSCRIPT)}"


def main() -> None:
    rows = list(csv.DictReader(CENSUS.open(), delimiter="\t"))
    if not rows:
        sys.exit(f"{CENSUS} is empty - run 22_closed_genome_search.py first")
    lengths = seq_lengths(FASTA)

    out = []
    for r in rows:
        acc = r["protein"]
        if acc not in lengths:
            sys.exit(f"{acc} is in the census but not in {FASTA.name}")
        out.append({
            "organism":       r["host"],
            "assembly":       r["assembly"],
            "replicon":       r["replicon"],
            "plasmid_name":   r["plasmid"],
            "nucleotide_acc": r["nucleotide"],
            "protein_acc":    acc,
            "length_aa":      lengths[acc],
            "pct_identity":   f"{float(r['pct_identity']):.1f} %",
            "evalue":         fmt_evalue(r["evalue"]),
            "_pid":           float(r["pct_identity"]),
        })

    # Plasmid block first, then chromosomal; each by descending identity.
    out.sort(key=lambda r: (r["replicon"] != "plasmid", -r["_pid"]))
    for i, r in enumerate(out, 1):
        r["rank"] = i
        del r["_pid"]

    with TABLE.open("w", newline="") as f:
        w = csv.DictWriter(f, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    n_plasmid = sum(r["replicon"] == "plasmid" for r in out)
    print(f"wrote {TABLE.relative_to(ROOT)}  "
          f"({len(out)} homologues: {n_plasmid} plasmid, "
          f"{len(out) - n_plasmid} chromosomal)")


if __name__ == "__main__":
    main()
