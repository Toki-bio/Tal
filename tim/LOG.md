# tim (Timema cristinae) — analysis log

*Scoped to this species only, per the same convention as `eri/LOG.md` and
`sco/LOG.md` — this is a separate, unrelated de novo SINE discovery
project sharing the SINEderella/report.html pipeline infrastructure with
the Talpidae species on this site, not a taxonomic relative. See
[`index.html`](../index.html)'s explicit "unconfirmed candidate" framing
for this section.*

---

## Why Timema, and why de novo

*Timema cristinae* (GCA_050494535.1, ASM5049453v1) is a stick insect
(Phasmatodea) — no known SINE consensus library exists for it, same
situation as the scorpion work above. Used the same de novo discovery
pipeline: `sine_scan.sh` (cross-species fragment homology search against
the 958-fragment SINEBase-derived query bank already prepared for the
scorpion project, reused as-is) to find candidate loci, followed by
manual `SubFam`-cluster-based consensus construction.

## De novo scan

Genome: 13 contigs, 1,227,621,598 bp — lightly fragmented, scanned fast
(unlike the worst scorpion genomes). Searched against the same
958-fragment query bank used for scorpions.

**Result: 924/958 query fragments (96.5%) got hits** — a substantially
higher hit rate than any scorpion genome — merging to **8,021** candidate
loci.

## Manual SubFam clustering (user's own work)

Ran `SubFam` with default parameters directly on the 8,021-sequence
candidate set, producing 161 chunk-consensus sequences. Tried the
automated `cluster_assist.js` first-pass tool as a shortcut: it found only
2 tiny clusters (9/161 sequences), 152 left as noise — matching its own
documented caveat that it does not reliably resolve real multi-subfamily
structure. Not usable; the user manually clustered the 161
chunk-consensuses into **8 candidate subfamily groups** (`t1`–`t8`, 2–9
underlying chunk-consensuses each).

## SINEderella classification run

`SINEderella GCA_050494535.1.fa candidates.q` (`run_20260821_132226`, full
5-step pipeline):

- Input: candidate loci from a fresh Step 1 search against the 8-sequence
  candidate consensus bank (independent from the raw 8,021-hit de novo
  scan pool the consensuses were clustered from).
- All 8 candidates received real genomic hits — copy counts ranged from
  t5 (1,010 copies) to t1 (172,236 copies).
- Full per-subfamily counts and stats in [`report.html`](report.html)'s
  Step2 section.

## Boundary refinement

Ran `step7_boundary_refine.sh` (the fraction-of-pairs-above-threshold
test, same method used for eri's e2-3 and all of scorpion) against this
run: **all 32 sides (8 subfamilies x 2 populations x 2 sides) confirmed**,
0 undetermined — a cleaner result than either eri or scorpion, likely
reflecting how distinct/well-separated the 8 candidate families are.
Numbers in [`boundary_refine/boundary_refinement.tsv`](boundary_refine/boundary_refinement.tsv).

top100 (bitscore-ranked) and general/rand100 populations were tested
separately per the same fix applied to scorpion (a boundary confirmed on
a random sample is not automatically valid for the bitscore-ranked set) —
see `report.html`'s two Extension columns per subfamily.

## Alignments

`step8a_extract_alignments.sh` + `step8b_publish_report.sh` ran cleanly
against this run with **zero warnings/failures** — all 8 subfamilies got
all three variants (top100/rand100/subfam), unlike scorpion where 6
subfamilies hit a real `bedtools getfasta` bug on certain loci. Consensus
sequences are first in every alignment (per the site-wide fix applied this
session).

**A real infrastructure bug found and worked around**: `SINEderella`'s own
orchestrator recreates `results/` as a symlink farm to other pipeline
output *after* `step6_report.sh` has already written a real
`results/report.html` file there, silently destroying it (confirmed
directly: the log showed "Wrote .../results/report.html (1.84 MB)" then
"Creating results/ directory" then completion, but the file was gone
afterward). Same class of issue previously worked around for scorpion by
supplying a locally-held copy; here, fixed by simply re-running
`step6_report.sh <run_root> --tal-species-code tim --no-sineplot`
standalone after the main orchestrator finished, which regenerates
`report.html` in the now-stable `results/` directory. Not yet fixed in
`SINEderella` itself — flagged as follow-up.

## What's here

- [`report.html`](report.html) — full report (Step1/Step2 stats,
  subfamily composition, bitscore thresholds, per-subfamily diagnostic
  plots, mutation-landscape PCA)
- `alignments/` — top100/rand100/subfam per subfamily, all with real
  genomic flanks (extension sizes vary — see report.html)
- `boundary_refine/boundary_refinement.tsv` — full dual-population
  boundary-refinement numbers

## 2026-08-21 — Rerun with refined consensuses (t1, t2, t345, t6, t7, t8)

User refined the 8 original candidate consensuses down to 6, based on
inspecting the first run's results — t3/t4/t5 merged into a single `t345`
consensus. Full `SINEderella` rerun (`run_20260821_145119`) against the
same genome with this refined 6-sequence bank:

| Subfamily | Assigned copies |
|---|---|
| t1 | 166,261 |
| t2 | 14,640 |
| t345 | 80,395 |
| t6 | 9,677 |
| t7 | 7,163 |
| t8 | 19,334 |

`step7_boundary_refine.sh`: all 24 sides (6 subfamilies x 2 populations x
2 sides) confirmed, 0 undetermined — same clean result as the first run.
`step8a`/`step8b`: 0 warnings, all 6 subfamilies got all three alignment
variants. Same `results/report.html`-gets-destroyed orchestrator bug hit
again (see below) — same workaround applied (standalone `step6_report.sh`
re-run).

This entirely replaces the previous 8-subfamily (`t1`-`t8`) publish —
old `t3`/`t4`/`t5` alignment files removed, all others regenerated fresh
against the new run.

## Not yet done

- Fix the `SINEderella` orchestrator's results/-directory race so
  `step6_report.sh` doesn't need a manual standalone re-run (hit again on
  this rerun, same as the first run — genuinely reproducible, not a fluke).
- Tandem-repeat contamination check (TRF) not done on this genome yet.
- No cross-referencing of the 6 refined candidate subfamilies against any
  Phasmatodea-specific literature (none was assumed to exist, same as
  scorpions — these remain unconfirmed candidates).
