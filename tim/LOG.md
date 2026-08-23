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

## 2026-08-23 — Cross-check against an independent de novo tool (AnnoSINE_v2)

A separate, unrelated exploration this session built a from-scratch
ab-initio SINE discovery pipeline and benchmarked several published tools
against it, including [AnnoSINE_v2](https://github.com/liaoherui/AnnoSINE_v2)
(Liao, Sun & Ou, *Mobile DNA* 2024) — a structure-based (TSD/poly-A/RNA-pol-III-promoter
signal) de novo SINE finder that needs no prior consensus library. Since it
needs no library, it's a genuinely independent check: it was never given the
6 `t1`/`t2`/`t345`/`t6`/`t7`/`t8` consensuses above, and found its own 55
candidates from the raw genome alone (mode 3 hybrid; zero hits from its
built-in homology/HMM library — Timema has no related family in it — so all
55 rest entirely on the structural module, each independently tagged as
tRNA-derived).

**Cross-check method**: `blastn` both directions between the 6 manually-curated
consensuses and AnnoSINE_v2's 55 candidates (`-evalue 1e-10`, `-word_size 11`).

**Result — real, mutual convergence, not circular**:

- **All 6 of the manually-curated subfamilies** have at least one significant
  AnnoSINE_v2 hit (e-values 1e-12 to 1e-95, 82–89% identity).
- **40 of AnnoSINE_v2's 55 candidates (73%)** independently converge onto one
  of the 6 subfamilies:

  | Subfamily | # AnnoSINE_v2 candidates matching it |
  |---|---|
  | t1 | 14 |
  | t2 | 17 |
  | t345 | 3 |
  | t6 | 6 |
  | t7 | 1 |
  | t8 | 2 |

  Identity is consistently 76–89% — high enough to clearly be the same
  underlying element, low enough (never >95%) that this reads as two
  independent measurements agreeing, not one derived from the other.
- **15 candidates (27%) found no match**: `SINE_3, 6, 17, 19, 22, 25, 28, 30,
  38, 39, 46, 47, 50, 52, 53`. Either genuinely novel families the manual
  6-subfamily curation hasn't captured, or false positives from AnnoSINE_v2's
  more permissive purely-structural net (no homology backstop at all for this
  taxon) — not yet manually inspected either way.

This is independent support that the 6 manually-curated subfamilies are real,
structurally-detectable SINEs and not artifacts of the `sear`/`ssearch36`
homology pipeline alone — a second, methodologically unrelated tool finds the
same families from raw sequence with no prior knowledge of them.

Full pairwise tables (on DRAGEN, not yet copied into this repo):
`/staging/tmp/timema_sines/cross_check/{user_vs_annosine,annosine_vs_user}.tsv`,
`unmatched.txt`. AnnoSINE_v2's raw 55-candidate output:
`/staging/tmp/timema_sines/annosine2/out/Seed_SINE.fa`.

## 2026-08-23 — `t7` set aside for separate treatment

`t7` is excluded from the new, further-split curated consensus set (below) —
not dropped as invalid, but held out deliberately. It only *partly* resembles
a SINE: the SINE-like signature appears to be a fragment embedded within a
larger repeat of unknown nature, not a clean standalone element the way
`t1`/`t2`/`t3`/`t6`/`t8` are. Classifying it alongside the others in a normal
SINEderella run would conflate two different questions (is this a SINE
copy-number/subfamily assignment vs. is this a fragment of a larger,
unrelated repeat family), so it needs its own separate investigation rather
than being folded into this run. Revisit once the surrounding larger repeat
is characterized.

## 2026-08-23 — Further manual split into finer subfamilies

User further split the curated consensus set into smaller, more homogeneous
subfamilies based on closer inspection — `t1` &rarr; `t1_1`/`t1-2`/`t1-3`/`t1-4`,
`t2` &rarr; `t2-1`/`t2-2`, `t345` split back apart into `t3-1`/`t3-2` (no longer
merged), `t6` &rarr; `t6-1`...`t6-5`, `t8` &rarr; `t8-1`/`t8-2` — **15 consensus
sequences total**, `t7` excluded (see above). Full `SINEderella` rerun
against this set (`run_20260823_063444`) completed successfully:

| Subfamily | Total assigned | Firm % | Sim mean |
|---|---|---|---|
| t1-2 | 120,953 | 28.97 | 0.6519 |
| t1_1 | 54,071 | 12.77 | 0.8135 |
| t3-1 | 47,215 | 12.20 | 0.9084 |
| t3-2 | 31,926 | 8.08 | 0.9091 |
| t1-4 | 40,902 | 7.95 | 0.3424 |
| t1-3 | 30,650 | 4.86 | 0.6126 |
| t8-2 | 21,133 | 4.20 | 0.6186 |
| t2-1 | 14,637 | 3.77 | 0.6852 |
| t8-1 | 11,429 | 2.83 | 0.6999 |
| t6-1 | 5,715 | 1.45 | 0.5626 |
| t6-2 | 2,720 | 0.65 | 0.6395 |
| t6-4 | 1,444 | 0.28 | 0.5585 |
| t6-3 | 883 | 0.16 | 0.5789 |
| t6-5 | 667 | 0.11 | 0.5840 |
| t2-2 | 381 | 0.06 | 0.8549 |

`t1-4` and several `t6-*`/`t2-2` subfamilies show noticeably lower `sim_mean`
and higher `leak_pct`/`conf_alt_pct` than the earlier 6/8-subfamily runs —
expected at this finer split (smaller, more homogeneous groups mean less
internal signal to firmly separate close neighbors, and boundary refinement
was not rerun at this granularity). Full numbers in
`results/summary.by_subfam.tsv` on DRAGEN and in the republished
[`report.html`](report.html).

`step8a_extract_alignments.sh` + `step8b_publish_report.sh` ran cleanly
against this run — all 15 subfamilies got top100/rand100 alignment variants,
14/15 also got the `subfam` variant (`t2-2` has only 215 assigned copies,
below the `MIN_COPIES_SUBFAM` threshold for SubFam re-clustering). This
**entirely replaces** the previous 6-subfamily (`t1`/`t2`/`t345`/`t6`/`t7`/`t8`)
publish — old alignment files removed, `tim_consensuses.fa` refreshed to the
new 15-sequence set, `report.html` regenerated fresh against this run. The
older `subfam/` and `boundary_refine/` directories are left in place as
historical record of the two earlier (8-subfamily and 6-subfamily) runs —
this v3 run did not rerun SubFam re-clustering site-wide or boundary
refinement, so there is nothing new to put there yet.

This was also the run that finally exercised the fixed `SINEderella`
orchestrator: `results/report.html` still got destroyed by the known
symlink-farm race a third time (worked around the same way as before,
standalone `step6_report.sh --no-sineplot` rerun), which is what motivated
actually fixing the race at the source afterward (block reorder in the main
`SINEderella` script, `results/` symlink farm now runs before Step 6, not
after) — not yet re-verified with a fresh run since the fix, but the fix
itself is in place for the next one.

## 2026-08-24 — t6-2/t6-3/t6-4 consensus extension, t2-1/t2-2 remerge (v4)

Two follow-ups from inspecting the v3 alignments directly in the MSA viewer:

**t6-2/t6-3/t6-4 consensuses were truncated.** The published alignments showed
a clearly shared, non-decaying motif continuing well past where the consensus
sequence stopped — not gradual decay into unique per-copy tails, a real
conserved block the original consensus simply never captured. Diagnosed by
column-level coverage analysis on the `rand100` alignment (excluding the old
consensus row): sustained shared signal continued for another 140–270
alignment columns past the old edge before genuinely dropping to background.
Fixed by re-extracting with a widened right flank (320bp, up from the
default 70bp) and re-aligning — not by rebuilding the consensus itself from
that alignment (a first attempt at that produced 40–60% `n` positions from
per-column vote splitting across a very indel-rich alignment, not usable).
User then manually extended the three consensuses directly using the wider
alignments.

**t2-1/t2-2 merged back into one `t2` consensus.**

Full `SINEderella` rerun against the updated 14-consensus set
(`run_20260823_103133`, "v4"):

| Subfamily | Total assigned | Firm % | Sim mean | vs. v3 |
|---|---|---|---|---|
| t1-2 | 120,871 | 29.02 | 0.6496 | ~same |
| t1_1 | 54,058 | 12.79 | 0.8830 | up from 0.8135 |
| t3-1 | 47,256 | 12.25 | 0.9208 | up from 0.9084 |
| t3-2 | 31,892 | 8.11 | 0.9201 | up from 0.9091 |
| t1-4 | 40,909 | 7.97 | 0.3452 | ~same |
| t1-3 | 30,666 | 4.89 | 0.6116 | ~same |
| t8-2 | 21,129 | 4.21 | 0.5793 | down from 0.6186 |
| t2 (merged) | 14,652 | 3.76 | 0.6488 | replaces t2-1/t2-2 |
| t8-1 | 11,400 | 2.83 | 0.6418 | down from 0.6999 |
| t6-1 | 5,705 | 1.45 | 0.5954 | ~same |
| t6-2 | 2,256 | 0.57 | 0.5589 | fewer copies (2,720→2,256), extended consensus is pickier |
| t6-4 | 1,061 | 0.27 | 0.5366 | fewer copies (1,444→1,061) |
| t6-3 | 717 | 0.17 | 0.5633 | fewer copies (883→717) |
| t6-5 | 957 | 0.10 | 0.6372 | ~same |

The extended t6-2/t6-3/t6-4 consensuses recruit fewer but more specific
copies, as expected from a longer, more precise reference. t1_1/t3-1/t3-2
all show a real `sim_mean` improvement from the refined consensus set.

`step8a`/`step8b` ran cleanly against this run — all 14 subfamilies got all
three alignment variants. This is the first run to complete after the
`SINEderella` orchestrator's results/-directory race fix (2026-08-23):
`results/report.html` (3.6MB) survived intact with no manual workaround
needed, confirming the fix holds.

This **entirely replaces** the v3 15-subfamily publish — `t2-1`/`t2-2`
alignment files removed, `t6-2`/`t6-3`/`t6-4` regenerated with the extended
consensuses, `tim_consensuses.fa` refreshed to the 14-sequence v4 set,
`report.html` regenerated fresh against this run.

## Not yet done

- Fix the `SINEderella` orchestrator's results/-directory race so
  `step6_report.sh` doesn't need a manual standalone re-run (hit again on
  this rerun, same as the first run — genuinely reproducible, not a fluke).
  **Fixed 2026-08-23, confirmed holding on the v4 run above.**
- Tandem-repeat contamination check (TRF) not done on this genome yet.
- No cross-referencing of the 6 refined candidate subfamilies against any
  Phasmatodea-specific literature (none was assumed to exist, same as
  scorpions — these remain unconfirmed candidates).
- Manually inspect the 15 AnnoSINE_v2 candidates with no match to any of the
  6 curated subfamilies (real/novel vs. structural-scan false positive).
- Characterize the larger repeat `t7` is embedded in, then decide how to
  treat it (separate family? chimera flag? exclude entirely?).
