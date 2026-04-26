# Tal SINE Analysis in Talpidae — Workflow

## Overview

Comprehensive analysis of Tal SINEs (tRNA-derived Short Interspersed Nuclear Elements)
in Talpidae (moles, desmans, shrew-moles). The workflow covers SINE detection, subfamily
classification, structural analysis (poly-A tails, terminators), and inter-species
comparison of orthologous insertion loci.

---

## Species

| Code | Species | Common name | Assembly |
|------|---------|-------------|----------|
| saq  | *Scalopus aquaticus* | Eastern mole | GCA_004024925.1 |
| ccr  | *Condylura cristata* | Star-nosed mole | GCF_000260355.1 |
| toc  | *Talpa occidentalis* | Western mole | TBD |
| gpy  | *Galemys pyrenaicus* | Pyrenean desman | GCA_019455555.1 |
| dmo  | *Desmana moschata* | Russian desman | TBD |

---

## Repository Layout

```
db/                    # Tal SINE consensus databases (.bnk, FASTA)
saq/                   # Scalopus aquaticus results
ccr/                   # Condylura cristata results
toc/                   # Talpa occidentalis results + alignments
gpy/                   # Galemys pyrenaicus results + alignments
dmo/                   # Desmana moschata results (pending)
scripts/               # Analysis scripts (Steps 5-8)
```

---

## Tools

- **SINEderella** — SINE search + subfamily assignment pipeline  
  Kit path: `/data/V/toki/bin/SINEderella/`  
  GitHub: https://github.com/Toki-bio/SINEderella
- **sear** — alignment-based SINE search (inside SINEderella)
- **asSINEment** — subfamily assignment (inside SINEderella)
- **MUSCLE / MAFFT** — multiple sequence alignment
- **minimap2 / LASTZ** — genome-genome alignment for liftover
- **liftOver** (UCSC) — coordinate mapping between assemblies

---

## Step 1 — SINE Search & Extraction

Uses `sear` (BLASTn-based) to search a query consensus bank against the genome.

**Input:** genome FASTA, `tal.bnk` (Tal consensus bank)  
**Output:** `extracted.fasta` (all hits ≥ threshold), `all_hits.labeled.bed`

```bash
# Kit setup (one species at a time)
BASE=/data/W/toki/Genomes/Mammalia/Eulipotyphla
SPECIES=saq   # ccr | toc | gpy | dmo

cd $BASE/$SPECIES
bash /data/V/toki/bin/SINEderella/step1_search_extract.sh \
    genome.fna tal.bnk 48
```

**Key parameters:**
- Chunk size: 30,000 bp
- Flank: 50 bp each side
- Threads: 48

---

## Step 2 — Subfamily Assignment (asSINEment)

All extracted copies are assigned to the best-matching subfamily consensus using
bitscores. Assignments are labeled firm/soft/unassigned.

**Input:** `extracted.fasta`, `tal.bnk`  
**Output:** `assignment_full.tsv`, `assignment_stats.tsv`, `summary.by_subfam.tsv`

```bash
bash /data/V/toki/bin/SINEderella/step2_asSINEment.sh
```

**Output columns (summary.by_subfam.tsv):**
`subfam | firm_assigned | soft_assigned | total_assigned | leak_n | conf_alt_n | firm_pct | total_pct | sim_mean | sim_median`

---

## Step 3 — Post-processing & Stats

Computes per-subfamily statistics, similarity score distributions,
self-bitscore normalisation.

```bash
bash /data/V/toki/bin/SINEderella/step3_postprocess.sh
```

**Key outputs:**
- `sim_scores.tsv` — pairwise similarity scores
- `self_bits_real.tsv` — normalised self-bitscores
- `regions.by_subfam.bed` — genomic coordinates by subfamily

---

## Step 4 — Plots & HTML Report

Generates interactive Plotly-based HTML report with:
- Subfamily size pie / bar charts
- Assignment confidence heatmaps
- Similarity score distributions
- Genomic density plots

```bash
bash /data/V/toki/bin/SINEderella/step4_plots.sh
# Output: results/report.html
```

---

## Step 5 — Tribe Identification (nearest-copy clusters)  *(TODO)*

Identify groups of near-identical copies that likely represent recent transposition bursts ("tribes").

**Approach:**
1. All-vs-all BLASTn within each subfamily (`assigned.fasta` per subfamily)
2. Retain pairs with identity ≥ 99%, alignment ≥ 90% of consensus length
3. Build similarity graph → MCL clustering (inflation I=2)
4. Output tribe table: tribe_id, copy_count, mean_identity, genomic_coords

```bash
# scripts/step5_tribes.sh  (to be implemented)
```

**Output:** `tribes.tsv` per species; columns: `copy_id | tribe_id | tribe_size | mean_identity`

---

## Step 6 — Poly-A Tail Analysis  *(TODO)*

Tal SINEs are transcribed by RNA Pol III and expected to have poly-A tails.

**Approach:**
1. Extract 50 bp 3′-downstream of each forward-strand SINE hit (BEDTools slop + getfasta)
2. Score poly-A: longest uninterrupted A-run, %A in window
3. Filter by strand (only + strand hits produce standard poly-A on reference)
4. Report fraction with poly-A ≥ 6 nt, length distribution per subfamily

```bash
# scripts/step6_polya.sh  (to be implemented)
```

**Output:** `polya_stats.tsv`; columns: `copy_id | strand | polyA_len | pctA_50bp | has_polyA`

---

## Step 7 — Terminator Analysis  *(TODO)*

RNA Pol III transcription terminates at a run of ≥ 4 T residues on the non-template strand.
Truncated copies may lack the terminator, leading to read-through transcription.

**Approach:**
1. Extract 20 bp immediately 3′ of each copy
2. Search for TTTT (or TTTTTT) motif within first 10 bp downstream
3. Also check for internal terminator (within copy body) in truncated copies
4. Cross-reference with copy length to classify full-length vs truncated

```bash
# scripts/step7_terminator.sh  (to be implemented)
```

**Output:** `terminator_stats.tsv`; columns: `copy_id | copy_len | has_3prime_term | internal_term_pos`

---

## Step 8 — Inter-species Comparison of Orthologous Loci  *(TODO)*

Identify Tal insertions that are shared (orthologous) vs. species-specific.

### 8a — Genome-genome alignment

```bash
# Reference species: saq (or toc when toc assembly is available chromosome-level)
minimap2 -a --cs -x asm5 \
    $BASE/saq/genome.fna \
    $BASE/ccr/genome.fna \
    > ccr_vs_saq.paf

# Convert PAF → chain/net for liftOver
```

### 8b — Liftover of SINE coordinates

```bash
# Map saq SINE BED → ccr coordinates
liftOver \
    $BASE/saq/results/regions.by_subfam.bed \
    saq_to_ccr.chain \
    saq_SINEs_in_ccr_coords.bed \
    unmapped.bed
```

### 8c — Orthologous locus table

For each successfully lifted-over locus:
- Retrieve ccr subfamily assignment at that coordinate
- Compare subfamily labels (conserved / switched / absent)
- Mark species-specific insertions (failed liftover = potentially lineage-specific)

```bash
# scripts/step8_orthologs.sh  (to be implemented)
```

**Output:** `ortholog_table.tsv`; columns:  
`locus_id | saq_subfam | ccr_subfam | toc_subfam | gpy_subfam | dmo_subfam | n_species | notes`

---

## Analysis Status

| Species | Step 1 | Step 2 | Step 3 | Step 4 | Step 5 | Step 6 | Step 7 | Step 8 |
|---------|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
| saq     | ✓ | ✓ | ✓ | ✓ | — | — | — | — |
| ccr     | ✓ | ✓ | ✓ | ✓ | — | — | — | — |
| toc     | ✓ | partial | — | — | — | — | — | — |
| gpy     | partial | — | — | — | — | — | — | — |
| dmo     | — | — | — | — | — | — | — | — |

---

## Kit Paths Reference

```
Genomes base:  /data/W/toki/Genomes/Mammalia/Eulipotyphla/
SINEderella:   /data/V/toki/bin/SINEderella/
saq genome:    .../saq/GCA_004024925.1_ScaAqu_v1_BIUU_genomic.fna
ccr genome:    .../ccr/GCF_000260355.1_ConCri1.0_genomic.fna
saq run:       .../saq/run_20260425_182219/
ccr run:       .../ccr/run_20260425_100809/
```

---

## Notes

- Tal SINEs are tRNA-derived; consensus length ~250–280 bp
- 9 subfamilies identified in saq (s1–s9), 7 in ccr (g2–g6)
- Subfamily naming convention follows SINEderella output (s = Scalopus; g = generic/ccr)
- For inter-species comparison, subfamilies need remapping to a shared nomenclature
- GWHHASN00000000.1 genome (at Eulipotyphla root) — species TBD, chromosome-level assembly
