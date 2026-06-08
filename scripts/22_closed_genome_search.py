#!/usr/bin/env python3
"""Closed-genome PL25_M23 homologue survey.

Pipeline:
  1. Parse the 250-hit BLASTp output (results/05_blast_nr/PL25_M23_vs_nr.tsv).
  2. For every halotolerant-Bacillaceae protein hit, query NCBI IPG to find
     source assemblies; keep only assembly_level ∈ {Complete Genome, Chromosome}.
  3. For every retained assembly, look up the source nucleotide-record title
     and assign replicon type (plasmid vs chromosome).
  4. Fetch full-length protein sequences for the retained homologues.

Outputs (all under results/22_closed_genome_search/):
  • closed_M23_homologues.tsv         — table of 14 closed-genome hits
  • closed_M23_homologues.faa         — 14 full-length protein sequences
  • closed_M23_with_outgroup.faa      — 16-sequence input for phylogeny
                                        (14 + PL25_M23 + lysostaphin)

Reproducibility:
  Run:  /usr/local/Caskroom/miniforge/base/envs/macsy/bin/python \\
          scripts/22_closed_genome_search.py
  Requires: Python ≥ 3.10, no third-party packages (urllib + xml only).
  Network: requires NCBI E-utilities access (~5 min for ~50 API calls).
"""
from __future__ import annotations
import sys, time, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

EMAIL  = "ayman.farran@kaust.edu.sa"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ROOT   = Path(__file__).resolve().parent.parent
BLAST  = ROOT / "results/05_blast_nr/PL25_M23_vs_nr.tsv"
OUTDIR = ROOT / "results/22_closed_genome_search"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Halotolerant Bacillaceae genera screened (15 total)
HALO_GENERA = {
    "Halobacillus","Lentibacillus","Virgibacillus","Gracilibacillus","Oceanobacillus",
    "Cytobacillus","Ornithinibacillus","Salinicoccus","Salimicrobium","Salibacillus",
    "Salinibacillus","Pontibacillus","Aquibacillus","Halolactibacillus","Salirhabdus",
}

# ────────────────────────────────────────────────────────────────────────────
def parse_blast(path: Path):
    """Yield (protein_acc, pct_id, evalue, organism, genus) per hit."""
    import re
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 14: continue
        acc = f[1].split("|")[1] if "|" in f[1] else f[1]
        pct, ev = f[2], f[10]
        m = re.search(r"\[([^]]+)\]", f[13])
        org = m.group(1) if m else ""
        genus = org.split()[0] if org else ""
        yield acc, pct, ev, org, genus

def fetch_ipg(acc: str):
    """Return [(assembly_acc, organism)] for a protein accession."""
    url = f"{EUTILS}/efetch.fcgi?db=protein&rettype=ipg&retmode=xml&id={acc}&email={EMAIL}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            root = ET.parse(r).getroot()
    except Exception:
        return []
    rows = []
    for cds in root.iter("CDS"):
        a = cds.attrib
        if a.get("assembly"):
            rows.append((a["assembly"], a.get("org",""), a.get("accver","")))
    return rows

def fetch_assembly_level(asm_acc: str) -> str:
    url1 = f"{EUTILS}/esearch.fcgi?db=assembly&term={asm_acc}&email={EMAIL}"
    try:
        with urllib.request.urlopen(url1, timeout=30) as r:
            uids = [e.text for e in ET.parse(r).getroot().iter("Id")]
        if not uids: return "unknown"
        url2 = f"{EUTILS}/esummary.fcgi?db=assembly&id={uids[0]}&email={EMAIL}"
        with urllib.request.urlopen(url2, timeout=30) as r:
            for e in ET.parse(r).getroot().iter("AssemblyStatus"):
                return e.text or "unknown"
    except Exception:
        return "error"
    return "unknown"

def fetch_nuccore_title(nuc_acc: str) -> str:
    url = f"{EUTILS}/esummary.fcgi?db=nuccore&id={nuc_acc}&email={EMAIL}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            for item in ET.parse(r).getroot().iter("Item"):
                if item.attrib.get("Name") == "Title":
                    return item.text or ""
    except Exception:
        pass
    return ""

def fetch_protein_fasta(acc: str) -> str:
    url = f"{EUTILS}/efetch.fcgi?db=protein&id={acc}&rettype=fasta&retmode=text&email={EMAIL}"
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()

def classify_replicon(title: str) -> tuple[str, str]:
    """Return (replicon_type, plasmid_name)."""
    t = title.lower()
    if "plasmid" in t:
        # try to extract plasmid name token after the word "plasmid"
        toks = title.split()
        if "plasmid" in toks:
            i = toks.index("plasmid")
            pname = toks[i+1].rstrip(",;:") if i+1 < len(toks) else "unnamed"
        else:
            pname = "unnamed"
        return "plasmid", pname
    if "chromosome" in t or "complete genome" in t:
        return "chromosome", "—"
    return "unknown", "—"

# ────────────────────────────────────────────────────────────────────────────
def main() -> None:
    print(f"Reading {BLAST.name}...", file=sys.stderr)
    halo_hits = [(acc, pct, ev, org)
                 for acc, pct, ev, org, gen in parse_blast(BLAST)
                 if gen in HALO_GENERA]
    print(f"  {len(halo_hits)} halotolerant Bacillaceae hits to check", file=sys.stderr)

    closed = []           # list of dicts: organism, assembly, replicon, plasmid, nuc, protein, pct, ev
    seen   = set()
    for i, (acc, pct, ev, org) in enumerate(halo_hits, 1):
        print(f"  [{i}/{len(halo_hits)}] {acc}", file=sys.stderr, end="\r")
        for asm, asm_org, nuc in fetch_ipg(acc):
            if asm in seen: continue
            seen.add(asm)
            if not asm.startswith("GCF_"): continue          # RefSeq only
            if (asm_org.split()[0] if asm_org else "") not in HALO_GENERA: continue
            time.sleep(0.34)
            if fetch_assembly_level(asm) not in ("Complete Genome","Chromosome"):
                continue
            time.sleep(0.34)
            title = fetch_nuccore_title(nuc)
            time.sleep(0.34)
            rep, pname = classify_replicon(title)
            closed.append(dict(
                organism=asm_org, assembly=asm, replicon=rep, plasmid=pname,
                nucleotide=nuc, protein=acc, pct_identity=pct, evalue=ev,
            ))
            print(f"\n  ✓ {asm}  {rep:11s}  {asm_org}", file=sys.stderr)
        time.sleep(0.34)
    print(file=sys.stderr)

    # Write TSV
    tsv = OUTDIR / "closed_M23_homologues.tsv"
    cols = ["organism","assembly","replicon","plasmid","nucleotide","protein","pct_identity","evalue"]
    with tsv.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in closed:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"✓ {tsv}  ({len(closed)} entries)", file=sys.stderr)

    # Write FASTA (just the 14 homologues)
    faa = OUTDIR / "closed_M23_homologues.faa"
    with faa.open("w") as out:
        for r in closed:
            org = r["organism"].replace(" ", "_")
            ctx = ("plasmid_"+r["plasmid"]) if r["replicon"]=="plasmid" else "chromosome"
            seq = "".join(fetch_protein_fasta(r["protein"]).strip().split("\n")[1:])
            out.write(f">{org}|{ctx}|{r['protein']}\n{seq}\n")
            time.sleep(0.34)
    print(f"✓ {faa}", file=sys.stderr)

    # Add PL25_M23 + lysostaphin outgroup (for phylogeny)
    seedset = ROOT / "results/07_M23_phylogeny/M23_seedset_v3.faa"
    combined = OUTDIR / "closed_M23_with_outgroup.faa"
    with combined.open("w") as out:
        keep = False; current_label = None
        for L in seedset.read_text().splitlines():
            if L.startswith(">"):
                if L.startswith(">PL25_M23"):
                    out.write(">PL25_M23|plasmid_pPL25-M23|this_study\n"); keep = True
                elif L.startswith(">REF_Lysostaphin"):
                    out.write(">Lysostaphin_outgroup|chromosome|P10547\n"); keep = True
                else:
                    keep = False
            elif keep:
                out.write(L + "\n")
        # append the 14 hits
        out.write(faa.read_text())
    print(f"✓ {combined}", file=sys.stderr)

if __name__ == "__main__":
    main()
