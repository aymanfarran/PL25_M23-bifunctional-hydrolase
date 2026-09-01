#!/usr/bin/env python3
"""Export the wet-lab antibacterial assay workbook to tidy TSV.

Reads the bench workbook (five sheets, one per assay) and writes one TSV per
assay into results/20_activity/. Those TSVs are the canonical raw-data record
kept in this repository; 41_activity_panels.py plots from them, not from the
workbook, so the figures can be rebuilt without the Excel file.

Assay conditions (all five): PL25_M23 vs OD600-normalised bacteria in PBS
pH 7.4 at 37 degC, 60 min unless the assay varies time; Miles-Misra drop-plate
enumeration; 3 technical replicates per condition, single biological experiment.

Usage:
    python scripts/40_activity_export_raw.py [path/to/data.xlsx]

Outputs:
    results/20_activity/turbidity_raw.tsv
    results/20_activity/dose_response_raw.tsv
    results/20_activity/time_kill_raw.tsv
    results/20_activity/edta_raw.tsv
    results/20_activity/salt_raw.tsv
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "20_activity"
OUT.mkdir(parents=True, exist_ok=True)

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / "Desktop" / "data.xlsx"

REPS = ["rep1_CFU_per_mL", "rep2_CFU_per_mL", "rep3_CFU_per_mL"]

# Workbook uses short strain labels; the TSVs carry the full binomials.
FULL_NAME = {
    "B. subtilis":      "Bacillus subtilis",
    "S. warneri":       "Staphylococcus warneri",
    "MRSA":             "MRSA",
    "V. salarius PL25": "Virgibacillus salarius PL25",
}
STRAINS = list(FULL_NAME)


def norm(v):
    """Strip non-breaking spaces the workbook uses inside strain labels."""
    return str(v).replace("\xa0", " ").strip() if pd.notna(v) else ""


def nums(row, cols):
    """Pull the three replicate CFU values out of a worksheet row."""
    return pd.to_numeric(row[list(cols)], errors="coerce").dropna().values[:3]


def fmt(vals):
    return [f"{v:.2e}" for v in vals]


def strain_rows(sheet, label_col):
    """Index the rows belonging to each strain, in worksheet order."""
    df = pd.read_excel(SRC, sheet, header=None)
    df[label_col] = df[label_col].map(norm)
    return df, {s: df.index[df[label_col] == s].tolist() for s in STRAINS}


# ── Turbidity / primary antibacterial assay (100 ug/mL, 60 min) ──────
def export_turbidity():
    df, idx = strain_rows("Turbidity ", 1)
    rows = []
    for s in STRAINS:
        untreated, treated = idx[s][0], idx[s][1]
        for label, i in (("Untreated", untreated), ("PL25_M23_100ugmL", treated)):
            rows.append([FULL_NAME[s], label, *fmt(nums(df.iloc[i], (3, 4, 5)))])
    return pd.DataFrame(rows, columns=["strain", "treatment", *REPS])


# ── Dose-response (0 / 50 / 100 / 300 ug/mL, 60 min) ─────────────────
def export_dose():
    DOSE_COLS = {0: (5, 6, 7), 50: (9, 10, 11), 100: (13, 14, 15), 300: (17, 18, 19)}
    df, idx = strain_rows("Dose killing", 3)
    rows = []
    for s in STRAINS:
        untreated, treated = idx[s][0], idx[s][1]
        for dose, cols in DOSE_COLS.items():
            for label, i in (("Untreated", untreated), ("PL25_M23", treated)):
                rows.append([FULL_NAME[s], label, dose, *fmt(nums(df.iloc[i], cols))])
    return pd.DataFrame(rows, columns=["strain", "treatment", "dose_ug_mL", *REPS])


# ── Time-kill (0-180 min at 100 ug/mL) ───────────────────────────────
def export_time_kill():
    TIME_COLS = {0: (4, 5, 6), 30: (8, 9, 10), 60: (12, 13, 14),
                 120: (16, 17, 18), 180: (20, 21, 22)}
    df, idx = strain_rows("Time killing", 2)
    rows = []
    for s in STRAINS:
        untreated, treated = idx[s][0], idx[s][1]
        for t, cols in TIME_COLS.items():
            for label, i in (("Untreated", untreated), ("PL25_M23", treated)):
                rows.append([FULL_NAME[s], label, t, *fmt(nums(df.iloc[i], cols))])
    return pd.DataFrame(rows, columns=["strain", "treatment", "time_min", *REPS])


# ── Metal dependence (1 mM EDTA +/- 100 uM ZnCl2, 100 ug/mL, 60 min) ─
def export_edta():
    CONDITIONS = {
        "Untreated":          "Untreated",
        "EDTA only":          "with 1mM EDTA only",
        "PL25_M23":           "PL25_M23 100 µg/mL",
        "PL25_M23+EDTA":      "PL25_M23 100 µg/mL + 1mM EDTA",
        "PL25_M23+EDTA+Zn":   "PL25_M23 100 µg/mL +1mM EDTA + ZnCl₂",
    }
    df = pd.read_excel(SRC, "EDTA + ZnCl2", header=None)
    df[2] = df[2].map(norm)
    df[3] = df[3].map(norm)
    rows = []
    for s in STRAINS:
        for out_label, sheet_label in CONDITIONS.items():
            hit = df[(df[2] == s) & (df[3] == sheet_label)]
            if hit.empty:
                continue
            rows.append([FULL_NAME[s], out_label, *fmt(nums(hit.iloc[0], (4, 5, 6)))])
    return pd.DataFrame(rows, columns=["strain", "condition", *REPS])


# ── Salt tolerance (0 / 0.5 / 1.0 M NaCl, 100 ug/mL, 60 min) ─────────
def export_salt():
    SALT = {
        0.0: "PL25_M23 100 µg/mL + 0 M NaCl",
        0.5: "PL25_M23 100 µg/mL + 0.5 M NaCl",
        1.0: "PL25_M23 100 µg/mL +1. 0 M NaCl",
    }
    df = pd.read_excel(SRC, "Salt ", header=None)
    df[3] = df[3].map(norm)
    df[4] = df[4].map(norm)
    rows = []
    for s in STRAINS:
        sub = df[df[3] == s]
        untreated = sub[sub[4] == "Untreated"]
        if untreated.empty:
            continue
        # The untreated control was measured once, at 0 M, and is the reference
        # for every salt condition of that strain.
        rows.append([FULL_NAME[s], "Untreated", 0.0,
                     *fmt(nums(untreated.iloc[0], (5, 6, 7)))])
        for nacl, sheet_label in SALT.items():
            hit = sub[sub[4] == sheet_label]
            if hit.empty:
                continue
            rows.append([FULL_NAME[s], "PL25_M23", nacl,
                         *fmt(nums(hit.iloc[0], (5, 6, 7)))])
    return pd.DataFrame(rows, columns=["strain", "treatment", "NaCl_M", *REPS])


EXPORTS = {
    "turbidity_raw.tsv":     export_turbidity,
    "dose_response_raw.tsv": export_dose,
    "time_kill_raw.tsv":     export_time_kill,
    "edta_raw.tsv":          export_edta,
    "salt_raw.tsv":          export_salt,
}

if __name__ == "__main__":
    if not SRC.exists():
        sys.exit(f"workbook not found: {SRC}")
    print(f"reading {SRC}\n")
    for name, fn in EXPORTS.items():
        df = fn()
        df.to_csv(OUT / name, sep="\t", index=False)
        print(f"  {name:<24} {len(df):>3} rows")
    print(f"\nwrote {len(EXPORTS)} files to {OUT}")
