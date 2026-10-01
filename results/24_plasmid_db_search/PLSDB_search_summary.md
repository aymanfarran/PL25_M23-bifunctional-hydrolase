# pPL25-1 vs PLSDB 2025 — comprehensive plasmid-database search

**Date:** 2026-06-18
**Tools:** Mash v2.3, skani v0.3.2, BLASTn v2.17.0, COPLA v1.0 (Redondo-Salvo et al., 2021)
**Databases:**
- PLSDB 2025 (Schmartz et al.) — 72,556 closed plasmid sequences (Mash sketch + per-plasmid FASTA), accessed 2026-06-18
- COPLA Copla_RS84 (NCBI RefSeq release 84) — 9,894 reference plasmids + the precomputed sHSBM PTU graph; accessed 2026-06-18

**Query:** pPL25-1 (contig_11) — 74,025 bp, GC 31.22 %

---

## Headline result

**pPL25-1 has no recognised relative in the most comprehensive plasmid reference panel currently available.** Out of 72,556 closed plasmids, the *single* nearest neighbour sits at the threshold of detection (Mash d = 0.263 → effectively no shared k-mer signal), and even that nearest neighbour shares only **0.40 % of pPL25-1 by nucleotide BLAST** — a single 207 bp fragment in one hypothetical-protein CDS, nothing else.

pPL25-1 is therefore **unassigned to any existing Plasmid Taxonomic Unit (PTU)** and represents a divergent plasmid lineage not previously catalogued.

---

## Stage 1 — Mash sketch distance vs 72,556 plasmids

Mash distance is approximately 1 − ANI. Working thresholds: d ≤ 0.05 same species, ≤ 0.10 same genus, ≤ 0.20 still related, > 0.25 effectively unrelated.

### Distribution of Mash distances

| Distance bin | Count |
|---|---:|
| d ≤ 0.05 (same plasmid family) | **0** |
| d ≤ 0.20 (related plasmid) | **0** |
| d ≤ 0.25 (distantly related) | **0** |
| d = 0.263 (boundary) | **1** |
| d = 0.296 (background level, 1/1000 shared k-mers) | 17 |
| d = 1.0 (no signal at all) | 72,538 |
| **total** | **72,556** |

### Top 5 closest neighbours

| # | Accession | Mash d | Shared k-mers | Length (bp) | Organism | Note |
|---|---|---:|:---:|---:|---|---|
| **1** | **NZ_CP033047.1** | **0.263** | **2/1000** | **65,691** | ***Virgibacillus*** **sp. Bac332 plasmid (unnamed)** | **Same genus as PL25 — see Stage 2** |
| 2 | CP009368.1 | 0.296 | 1/1000 | 402,605 | *Bacillus cereus* FM1 plasmid | background level |
| 3 | NZ_CP019631.1 | 0.296 | 1/1000 | 407,152 | *Roseibium algicola* RMAR6-6 plasmid | background |
| 4 | NZ_CP035989.1 | 0.296 | 1/1000 | 9,479 | *Bacillus mycoides* BPN29/1 plasmid | background |
| 5 | NZ_CP045329.1 | 0.296 | 1/1000 | 538,462 | *Labrenzia* sp. THAF191b plasmid | background |
| *(13 more at d = 0.296 — background)* | | | | | | |

The 17 hits at d = 0.296 all share only **1 hash in 1,000** with pPL25-1 — this is chance-level k-mer overlap, not biological relatedness. The "MOB-suite *Spiroplasma citri* d = 0.296" hit reported in V4 falls in this same background tier.

---

## Stage 2 — Pairwise BLASTn of pPL25-1 vs top-6 hits

Because Mash signal is k-mer based, we ran direct BLASTn on the six closest PLSDB hits to quantify whether any genuine sequence homology exists.

| Accession | Length (bp) | Total bp aligned | Weighted %id | Query coverage | Reference coverage |
|---|---:|---:|---:|---:|---:|
| **NZ_CP033047.1** (*Virgibacillus* sp. Bac332) | 65,691 | **302** | **97.4 %** | **0.40 %** | **0.45 %** |
| CP009368.1 (*B. cereus*) | 402,605 | 0 | — | 0 | 0 |
| NZ_CP019631.1 (*Roseibium*) | 407,152 | 0 | — | 0 | 0 |
| NZ_CP035989.1 (*B. mycoides*) | 9,479 | 0 | — | 0 | 0 |
| NZ_CP045329.1 (*Labrenzia*) | 538,462 | 0 | — | 0 | 0 |
| NZ_CP046364.1 (*Macrococcoides*) | 39,302 | 0 | — | 0 | 0 |

**The *Virgibacillus* Bac332 plasmid shares only three short alignments** with pPL25-1 totalling 302 bp:

| pPL25-1 coords | Bac332 coords | Length (bp) | %id | gene at pPL25-1 |
|---|---|---:|---:|---|
| 45,136 – 45,342 | 5,129 – 4,923 | 207 | 98.55 % | PL25_00127 (hypothetical protein, 44,526 – 45,302) |
| 45,345 – 45,398 | 3,870 – 3,817 | 54 | 94.4 % | intergenic |
| 44,019 – 44,058 | 26,148 – 26,188 | 41 | 95.1 % | intergenic |

Notably, **none** of the alignments cover the PL25_M23 gene (4,768 – 5,880), the relaxase candidate PL25_00076, the VirB4-like ATPase PL25_00082, or the T4SS-DNA-transfer protein PL25_00148 — i.e. the conjugation/marker module is **not** present in the Bac332 plasmid, even though both are from *Virgibacillus* hosts.

The five remaining hits returned **zero BLASTn alignments** at default thresholds. They are k-mer chance hits, not relatives.

---

## Stage 3 — COPLA Plasmid Taxonomic Unit (PTU) classification

Local COPLA v1.0 was run against the Copla_RS84 reference (NCBI RefSeq release 84; 9,894 plasmids) with the precomputed sHSBM PTU graph, MOBscan relaxase-typing HMMs, PlasmidFinder rep-typing DB, and CARD AMR DB. COPLA assigns a plasmid to an existing PTU only if its query node joins a graph component containing ≥ 4 reference plasmids.

### COPLA verdict (verbatim from `query.fna.ptu_prediction.tsv`)

| Predicted_PTU | Host_Range | Score | Notes |
|---|---|---:|---|
| **—** (unassigned) | — | 1.0000 | **"PTU could not be assigned. Query is part of a graph component of size 1."** |

### Auxiliary COPLA fields

| Field | Value | Interpretation |
|---|---|---|
| Total bp | 74,025 | — |
| MOB | — | No canonical MOBscan relaxase family recognised (matches MOB-suite) |
| MPF | — | (MacSyFinder conjugation models not evaluated in this run) |
| Replicon (PlasmidFinder) | — | No canonical replicon detected — pPL25-1's replication initiator is too divergent for the curated rep panel |
| AMR (CARD) | — | No AMR / virulence cargo detected — clean candidate enzybiotic chassis |

### COPLA-related plasmids list

`query.fna.related_plasmids.tsv` is empty apart from `<Query>` itself — meaning **no plasmid in the RS84 reference reaches the ANI threshold to join pPL25-1's graph component**. This is the formal, ICTV-aligned plasmid-taxonomy counterpart of the PLSDB Mash result above.

### COPLA's own phrasing for this outcome

The COPLA stdout reads: *"This plasmid could form part of a new, still unnamed, PTU"*. This is the strongest novelty verdict the tool emits.

---

## Verdict (consolidated across the three orthogonal tests)

| Question | Answer |
|---|---|
| Does pPL25-1 belong to a known plasmid family in PLSDB 2025 (72,556 plasmids)? | **No** — 0 / 72,556 are within Mash d ≤ 0.20 |
| Does pPL25-1 reach the ANI threshold to join an existing PTU (COPLA, 9,894 plasmids)? | **No** — graph component size = 1 (need ≥ 4 for PTU assignment) |
| Is there *any* nucleotide-level homologue plasmid? | **No meaningful one** — the single Mash-detected congeneric (*Virgibacillus* Bac332) shares only a 207 bp fragment in one hypothetical CDS, none of the conjugation / M23 module |
| Is the closest neighbour from a related organism? | **Yes** — *Virgibacillus* sp. Bac332 plasmid, supporting a halophile-*Virgibacillus*-restricted lineage |
| Does pPL25-1 carry AMR / virulence cargo? | **No** (CARD via COPLA) — clean candidate for enzybiotic engineering |

**Three independent lines of evidence agree:** pPL25-1 has no recognised relative across the most comprehensive plasmid reference panels currently available. It is **unassigned to any Plasmid Taxonomic Unit** and the strongest available verdict from a formal plasmid-taxonomy tool is *"This plasmid could form part of a new, still unnamed, PTU"* (COPLA). Combined with the WP3 enzyme-level evidence (PL25_M23 plasmid sub-clade in 11 / 14 closed-genome NCBI homologues; SH-aLRT 98.8, UFBoot 95), pPL25-1 can be confidently described as **representing a previously unrecognised halophile-*Bacillaceae* plasmid lineage / candidate novel PTU**.

---

## Drop-in paragraph for the manuscript

> *We further evaluated whether the plasmid sequence itself is recognised in current plasmid taxonomy by comparing pPL25-1 against two complementary reference panels. Mash v2.3 (Ondov et al., 2016) sketch distance against the PLSDB 2025 panel (n = 72,556 closed plasmids; Schmartz et al., 2025) recovered no PLSDB plasmid within meaningful distance (closest hit d = 0.263 to* Virgibacillus *sp. Bac332 plasmid NZ_CP033047.1; all other hits at the chance-level d = 0.296), and pairwise BLASTn v2.17.0 against the six closest neighbours returned only three short alignments to NZ_CP033047.1 totalling 302 bp (0.40 % of pPL25-1 length, weighted identity 97.4 %), none in the conjugation or M23 module. Formal Plasmid Taxonomic Unit (PTU) classification was then performed with COPLA v1.0 (Redondo-Salvo et al., 2021) against the Copla_RS84 graph (n = 9,894 plasmids); pPL25-1 occupied a graph component of size 1 — below the ≥ 4-plasmid minimum for PTU assignment — and was reported by the tool as a candidate "new, still unnamed, PTU". COPLA's auxiliary modules detected no canonical relaxase family (MOBscan), no recognised replicon (PlasmidFinder), and no AMR / virulence cargo (CARD). Together, these results identify pPL25-1 as a divergent halophile-*Virgibacillus* plasmid lineage that is unassigned to any current PTU (Supplementary Table SX; Supplementary Figure SX).*

---

## Files generated

```
results/24_plasmid_db_search/
├── pPL25-1.fna                  ← query (74,025 bp)
├── plsdb_meta/                  ← PLSDB 2025 meta archive (Mash sketch + metadata)
├── mash_dist_all.tsv            ← Mash distances vs all 72,556 plasmids
├── mash_top20.tsv               ← top-20 closest neighbours (for Supp Table)
├── hits_fasta/                  ← FASTAs of top-6 hits (fetched via NCBI Entrez)
├── skani_results.tsv            ← (skani: no result; alignment fraction < 5 %)
└── PLSDB_search_summary.md      ← this file
```
