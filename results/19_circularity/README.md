# contig_11 circularity check — FASTA-only fallback

**Status: inconclusive without long-read FASTQ.**

## What was attempted

| Test | Result |
|---|---|
| Python self-overlap: first *k* bp vs last *k* bp (exact) for *k* ∈ {2000, 1500, 1000, 500, 200, 100, 50} | No exact overlap ≥ 50 bp |
| Python self-overlap: ≥ 95 % identity over ≥ 200 bp | Best ≈ 32 % identity (random baseline) |
| Slide 5' 200-mer across the last 5 kb (approximate match) | Best ≈ 40 % identity (essentially random) |
| `nucmer --maxmatch --nosimplify -l 20 -c 50` (MUMmer 4.0.1) self-alignment | Only the trivial self-self diagonal (100 % over full 74,025 bp); **no terminal duplication of any length** |

## Why this is inconclusive

Long-read assemblers, including Flye, typically **trim the small (~200–2000 bp) circular-closure overlap after assembly** so the final contig is the "deduplicated" circular sequence with neither end repeated. A trimmed circular plasmid and a genuinely linear contig therefore look identical at the FASTA level.

The definitive answer is in **`assembly_info.txt`** (the file Flye writes alongside `assembly.fasta`). That file has one line per contig with a `circular` field (Y / N) and is the gold-standard record.

## What to run when the original Nanopore reads can be located

Either (a) recover `assembly_info.txt` and read the `circular` flag for contig_11, **or** (b) re-run the closure test with the original reads using Circlator on IBEX:

```bash
# On IBEX, in a node with ≥ 8 CPU, ≥ 16 GB RAM
module load mamba   # or conda
mamba create -y -n circlator -c bioconda -c conda-forge circlator mummer4 bwa samtools spades canu
mamba activate circlator

# Inputs:
#   contig_11.fna            ← the 74 kb assembled contig
#   PL25_nanopore.fastq.gz   ← original Nanopore reads
mkdir -p circ_run
circlator all                              \
  --threads 8                              \
  --assembler canu                         \
  --data_type nanopore-raw                 \
  contig_11.fna                            \
  PL25_nanopore.fastq.gz                   \
  circ_run

# Key outputs to check:
#   circ_run/06.fixstart.fasta      — final closed-circle FASTA, if circularised
#   circ_run/04.merge.circularise.log
#   circ_run/05.clean.contigs.fa
```

If Circlator successfully circularises contig_11, the log will explicitly report:
`Circularising contig contig_11`  followed by  `Circularised`.

## Provenance of the files in this directory

- `contig11_self.delta` / `contig11_self.coords` — MUMmer nucmer self-alignment (24 May 2026)
- This README — written as a manuscript-provenance record
