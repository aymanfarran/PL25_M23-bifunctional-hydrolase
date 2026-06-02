#!/usr/bin/env python3
"""
oriT consensus-motif scan on PL25 contig_11.

Implements the exploratory motif screen described in the manuscript
Methods and Supplementary Methods S4. Five canonical, primary-reference-
supported oriT nick-site motifs are tested on both strands of the input
sequence by exact match. The script also reports ≥ 18 bp exact inverted
repeats within a configurable window of the predicted relaxase candidate.

The screen is exploratory; positive hits would require experimental
confirmation of nicking and transfer before being reported as an oriT.

Usage:
    python scripts/11_orit_motif_scan.py [fasta] [relaxase_centre_bp]

Defaults to PL25 contig_11 (74,025 bp) and a relaxase centre at 5,500 bp
(approximate centroid of the predicted replication-relaxation-family
protein PL25_00076).
"""
import sys
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# Default inputs
ROOT = Path("/Users/farrana/prophage-endolysin-pipeline/non_endolysin_paper")
DEFAULT_FASTA  = ROOT / "results/03_M23_ICE_locus/contig_11.fna"
DEFAULT_CENTRE = 5500   # bp position of PL25_00076 (replication-relaxation candidate)
IR_HALF_WINDOW = 3000   # ± bp around the relaxase centre to scan for inverted repeats
IR_MIN_LEN     = 18     # minimum inverted-repeat length (bp)
IR_MAX_SEP     = 200    # maximum spacing between the two arms of the inverted repeat

# Canonical oriT nick-site motifs — each with a primary reference (see
# Supplementary Methods S4 of the manuscript for citations).
MOTIFS = {
    "IncP_RP4_nick"             : "CTGCCTGGTGCT",   # Pansegrau & Lanka 1996, NAR 24:4691
    "IncQ_RSF1010_nick"         : "CCGCCTGGTACCC",  # Scherzinger et al. 1993, EJB 217:929
    "IncW_R388_nick"            : "CTAACCGTGCAC",   # Llosa et al. 1994, J Mol Biol 235:448
    "pIP501_GramPos_nick"       : "TAGTGCGCC",      # Kurenbach et al. 2003, Plasmid 50:86
    "ICEBs1_nic_core_BACSU"     : "TATCATTAATAAT",  # Lee, Babic & Grossman 2010, MolMicro 75:268
}


def scan_motifs(seq: str) -> dict:
    """Exact-match scan of each motif on the forward strand and the
    reverse complement; returns a dictionary of hits."""
    rc = str(Seq(seq).reverse_complement())
    hits = {}
    for name, motif in MOTIFS.items():
        hits[name] = {"fwd": [], "rev": []}
        # Forward strand
        for m in re.finditer(motif, seq):
            hits[name]["fwd"].append(m.start() + 1)  # 1-based
        # Reverse strand
        for m in re.finditer(motif, rc):
            fwd_pos = len(seq) - (m.start() + len(motif)) + 1   # 1-based on fwd
            hits[name]["rev"].append(fwd_pos)
    return hits


def scan_inverted_repeats(seq: str, centre: int, half_window: int,
                          min_len: int, max_sep: int):
    """Look for exact inverted repeats of length ≥ min_len within the
    window [centre-half_window, centre+half_window], with the two arms
    separated by ≤ max_sep bp."""
    lo, hi = max(0, centre - half_window), min(len(seq), centre + half_window)
    window = seq[lo:hi]
    found = []
    # Search progressively shorter IRs; stop at IR_MIN_LEN
    for L in range(40, min_len - 1, -2):
        for i in range(len(window) - L * 2):
            arm1 = window[i:i + L]
            if "N" in arm1:
                continue
            rc_arm1 = str(Seq(arm1).reverse_complement())
            j = window.find(rc_arm1, i + L, i + L + max_sep)
            if j >= 0:
                found.append({
                    "len_bp"      : L,
                    "arm1_pos"    : lo + i + 1,        # 1-based on full sequence
                    "arm2_pos"    : lo + j + 1,        # 1-based
                    "arm1_seq"    : arm1,
                    "spacer_bp"   : j - (i + L),
                })
                if len(found) >= 10:
                    return found
    return found


def main():
    fa     = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_FASTA
    centre = int(sys.argv[2])  if len(sys.argv) > 2 else DEFAULT_CENTRE
    rec    = next(SeqIO.parse(str(fa), "fasta"))
    seq    = str(rec.seq).upper()

    print(f"Sequence            : {rec.id}  ({len(seq):,} bp)")
    print(f"Relaxase centre (bp): {centre:,}  (IR scan window ±{IR_HALF_WINDOW:,} bp)")
    print(f"Motifs tested       : {len(MOTIFS)} (primary-reference-supported)")
    print()

    # -- Motif hits --
    hits = scan_motifs(seq)
    any_hit = False
    for name, h in hits.items():
        if h["fwd"] or h["rev"]:
            any_hit = True
            for pos in h["fwd"]:
                print(f"  HIT  {name:25s}  fwd  pos {pos:>7,}  motif {MOTIFS[name]}")
            for pos in h["rev"]:
                print(f"  HIT  {name:25s}  rev  pos {pos:>7,}  motif {MOTIFS[name]}")
    if not any_hit:
        print("  No canonical oriT motif detected on either strand.")
    print()

    # -- Inverted repeats near relaxase --
    print(f"Inverted-repeat search (≥ {IR_MIN_LEN} bp, ≤ {IR_MAX_SEP} bp spacer):")
    irs = scan_inverted_repeats(seq, centre, IR_HALF_WINDOW, IR_MIN_LEN, IR_MAX_SEP)
    if irs:
        for ir in irs:
            print(f"  IR  L={ir['len_bp']:>2}  pos {ir['arm1_pos']:>7,}–{ir['arm2_pos']:>7,}  "
                  f"spacer {ir['spacer_bp']:>3} bp  arm1 5'-{ir['arm1_seq']}-3'")
    else:
        print(f"  No inverted repeat of ≥ {IR_MIN_LEN} bp detected near the relaxase candidate.")

    print()
    print("NOTE: This screen is exploratory. Positive motif hits require experimental")
    print("confirmation of nicking and conjugative transfer before being reported as an oriT.")


if __name__ == "__main__":
    main()
