# sco (Scorpiones) — analysis log

*Scoped to this species only, per the same convention as `eri/LOG.md` —
this is a separate, unrelated de novo SINE discovery project sharing the
SINEderella/report.html pipeline infrastructure with the Talpidae species
in this repo, not a taxonomic relative. Full methodology, all 11 genomes'
timing/loci counts, and the fragmentation-bottleneck engineering work are
in the SINEderella project's own `SCORPION_WORKFLOW.md` — this file covers
just this genome's classification run and links to that source.*

---

## Why scorpions, and why de novo

No known SINE consensus library exists for Scorpiones, so `SINEderella`
alone can't be used directly (it only classifies against already-known
consensuses). Used a separate de novo discovery step first
(`sine_scan.sh`, cross-species fragment homology search against a
958-fragment SINEBase-derived query bank) to find candidate loci, per
genome, across all 11 Scorpiones assemblies on NCBI (8 species — explicit
scope decision: every assembly, every version, not deduped to
one-per-species).

## Genome-wide bottleneck and its fix

The original all-11-genomes batch got stuck 7+ hours on one 720,000-contig
genome while easier, often-larger genomes sat unscanned. Root cause:
fragmentation (contig count) drives `sine_scan.sh`'s per-contig
`samtools faidx` subprocess-spawn overhead, not genome size — confirmed
directly (the two largest genomes by file size have the fewest contigs and
scan fastest). Fixed two ways: reordering the batch easy-first by contig
count with an adaptive time-budgeted per-genome escalation, and (for the
two worst genomes) concatenating all scaffolds into one contiguous
sequence before scanning, eliminating the fragmentation bottleneck
entirely (validated: 6.6x speedup on a 10% test subset, ~9-10% chimeric
hit rate correctly flagged/excluded via a coordinate-offset remap). All 11
genomes scanned; full per-genome numbers in `SCORPION_WORKFLOW.md`.

## This genome: GCA_982267015.1 (Belisarius xambeui)

Picked for further work as the run with the most raw hits (22,335 merged
loci from the de novo scan — largest of the 11).

**Manual SubFam clustering** (user's own work, not automated): ran `SubFam`
with default parameters directly on the full 22,335-sequence hit set,
producing 447 chunk-consensus files. Manually clustered these into **29
candidate subfamily groups** (`g01`–`g29`), each group's size (2–14
underlying chunk-consensuses) recorded in the consensus names.

**SINEderella classification run** (`SINEderella GCA_982267015.1.fa
candidates.q`, `run_20260821_025634`, THREADS=64, CHUNK_BP=30000,
FLANK=50, MODE=full):

- Input: 251,164 extracted candidate loci (Step 1, full-genome de novo
  search — note this re-searched with the 29-consensus bank as query,
  giving a different/larger candidate pool than the raw 22,335-hit set
  SubFam was clustered from)
- Unanimous 10/10 votes: 249,088
- **Assigned (passed threshold): 169,262 (67.4%)**
- Unassigned: 81,902 (32.6%)

Per-subfamily assignment ranges from `g06` (52,909 copies, the single
largest) down to `g20` (72 copies) — full breakdown in
[`report.html`](report.html)'s Step2 section.

## A second real bug found and fixed before publishing

Direct verification (not trusting the extraction log) caught a real bug:
5/20 sampled records in one subfamily's alignment were unflanked despite
the log saying "OK ... 50L+70R flanked" — all 5 had a malformed strand
notation `(+,-)` in their header (a merged/overlapping-hit artifact, same
class seen in `eri`'s data). The header-parsing regex only matched a
single-character strand suffix (`\([+-]\)$`), so `(+,-)` was left
un-stripped — and since it contains a literal `-`, it corrupted the
downstream coordinate split, silently dropping that BED line (no error,
just missing from the flanked output, which then fell back to the
unflanked raw sequence with no visible sign anything had gone wrong).

Fixed by stripping ANY trailing parenthetical rather than requiring
exactly one character, in both `extract_top100_rand100_subfam.sh` and
`extract_tribes.sh`. Re-verified on the corrected output: 600 sampled
records across 6 files, only 1 mismatch — a genuine scaffold-edge
truncation (locus starts at position 1 of its contig, no room for a full
flank on one side), not a bug. 117 of the 600 checked records specifically
had `(+,-)` headers and all passed cleanly after the fix.

**This same bug most likely affects `eri`'s already-published
top100/rand100/subfam alignments too** (eri's `assigned.fasta` has
~15,390 `(+,-)`-headers out of 1.76M, a similar ~0.9% rate) — not yet
re-verified/re-published there. Flagged as follow-up work, not done in
this session.

## What's here

- [`report.html`](report.html) — full report (Step1/Step2 stats,
  subfamily composition, bitscore thresholds, quality flags, per-subfamily
  diagnostic plots, mutation-landscape PCA)
- `alignments/` — top100 (by bitscore) / rand100 (random sample) / subfam
  (SubFam re-clustering) alignments per subfamily, all with real
  50bp-left/70bp-right genomic flanks (same convention as the Eulipotyphla
  species pages) — generated with the same
  `extract_top100_rand100_subfam.sh` used for `eri`, adapted for DRAGEN
  (KIT-only tool paths made overridable via env vars; DRAGEN has no
  equivalent to KIT's `/data/V/toki/bin/sample` tool, so an inline
  shuf-based fallback was added — first version of that fallback had a
  real bug (`cat | awk | shuf | head | awk` under `set -o pipefail`: `head`
  closing early sends SIGPIPE back up the pipe, which `pipefail` treats as
  failure and `set -e` then kills the whole script silently, mid-run, no
  error logged) — fixed by writing intermediate files instead of piping
  directly into `head`, matching a pattern already documented elsewhere in
  this codebase that this script's own author had apparently seen before
  but missed applying to the new fallback code.

## Not yet done

- Cross-referencing whether the earlier tRNA-gene-contamination hypothesis
  (raised during the single-genome exploratory work, before this batch —
  see `SCORPION_WORKFLOW.md`) applies to any of these 29 subfamilies
  specifically — not tested against this data yet.
- No tandem-repeat contamination check (TRF, the same one run for `eri`)
  has been done on this genome's assignments yet.
- The other 10 scanned genomes have raw de novo hits but no SINEderella
  classification run against a consensus bank yet — this genome was
  picked first because it had the most raw hits, not because the others
  are less interesting.
