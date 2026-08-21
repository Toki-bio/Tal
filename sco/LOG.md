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

## 2026-08-21 — Re-extraction with the `(+,-)` fix, verified before publishing

Re-ran `extract_top100_rand100_subfam.sh` (with the strand-parsing fix
described above) against this run — 83 alignment files, ~4.5 min. Verified
directly before publishing, not by trusting the log: 39 records sampled
across 14 files, weighted toward records with `(+,-)` headers (27 of the 39
— the exact case that failed before), checked
`ungapped_length == core_length + 120` (50L+70R). **0 mismatches.** Also
scanned all 83 files for `@U@` leaks: **0 occurrences** anywhere. This
confirms the fix actually works, not just that it runs without error.

`subfam/` rebuilt from the 25 `*_subfam.aln.fa` outputs (matches the 25
subfamilies with >=400 members; g20/g24/g26/g28 correctly skipped).

**Same bug very likely still affects `eri`'s already-published
top100/rand100/subfam alignments** — eri's `assigned.fasta` has a similar
`(+,-)` header rate and was extracted before this fix existed. Not
re-verified or re-published there in this session; flagged again as
follow-up work, not started.

## 2026-08-21 — Ran the modular boundary-refinement step against all 29 subfamilies

First real test of `step7_boundary_refine.sh` (the standalone, independently-
runnable module underlying the eri e2-3 boundary-extension work above —
same fraction-of-pairs-above-threshold test, same 45% cutoff, calibrated
against this genome's own background) against a full multi-subfamily run
rather than a single hand-picked subfamily. Deployed to DRAGEN
(`/staging/tmp/step7_boundary_refine.sh`), run directly against
`run_20260821_025634` (needed `setfacl -R -m user:copilot:rwx` granted by
the user first — the run directory is root-owned).

**Real bug found on this first execution** (never caught by `bash -n`,
since it's an awk-level error, not a bash syntax error): the per-subfamily
member-filtering step used `awk -v sub="$subfam" '$1 == sub'`, and `sub`
collides with awk's own built-in `sub()` function name —
`awk: fatal: cannot use gawk builtin 'sub' as variable name`, killing the
script immediately after Step A. Fixed by renaming the variable to `sf`.
Fixed in the SINEderella worktree, not yet pushed/merged upstream.

**Result**: completed in 35 seconds for the whole genome (29 subfamilies x
2 sides = 58 rows). 55/58 sides confirmed (boundary found within the
1000bp cap), 3 undetermined (g17 downstream, g18 upstream, g22 upstream —
still above the elevated-fraction threshold even at 1000bp, meaning either
a genuinely long fuzzy boundary or tandem/repeat contamination extending
further than tested). Confirmed boundaries range from 50bp (many
subfamilies — base flank already unique) up to 1000bp (g21 downstream).
Full table: [`boundary_refine/boundary_refinement.tsv`](boundary_refine/boundary_refinement.tsv).

Confirms the tool works as a standalone module against a real multi-
subfamily run, not just the single hand-tested eri e2-3 case — this was
the point of the test (per the project's modular-pipeline design: each
SINEderella step should be runnable independently, not just as part of
the full orchestrator).

**Not yet done**: results have not been fed into step8a/8b to regenerate
extended alignments or update `report.html` for scorpion — this run was a
tool-correctness test, not yet a publish step. The 3 undetermined
subfamilies have not been individually inspected by eye.

## 2026-08-21 — Published: extended alignments for all 29 subfamilies

Ran `step8a_extract_alignments.sh` + `step8b_publish_report.sh` against
this run, using `boundary_refinement.tsv` above to size each subfamily's
flanks (base 50bp/70bp + the confirmed per-side extension) instead of the
fixed 50L/70R every other species page on this site still uses.

**A second real bug found on this first execution**: the `subfam` variant's
element-extraction `bedtools getfasta` call was not wrapped in `set +e/-e`
like every other bedtools/mafft call in this script — when it failed for
`g02` (30284 members, the second subfamily processed), `set -e` killed the
whole script silently, with the real error swallowed by `2>/dev/null` and
no WARNING logged (the temp-dir cleanup trap fired cleanly, making it look
like a graceful finish rather than a crash). Only caught by rerunning with
`bash -x` tracing to reproduce and pinpoint the exact failing command.
Fixed by wrapping the call in `set +e/-e`, capturing stderr instead of
discarding it, and skipping just that subfamily's `subfam` variant (with a
logged WARNING) instead of killing the whole run. 6 subfamilies (g02, g05,
g07, g09, g10, g13) hit this on the real data and now skip gracefully
rather than crash.

**A publishing-integration bug also found and fixed**: `report.html`
already had a real (not placeholder) "Subfamily Alignments" table from an
earlier manual publish, so `step8b`'s placeholder-search-and-replace logic
didn't find a match and instead appended a second, duplicate table before
`</body>` — leaving the page with two `id="alignments"` sections and the
old one's hardcoded "50bp upstream + 70bp downstream" text now factually
wrong (flanks vary per subfamily since boundary refinement). Fixed by hand:
removed the old section, restyled the new step8b-generated table to match
the site (card styling, intro paragraph, added an Extension column showing
each subfamily's confirmed up/down boundary).

**A verification false alarm, resolved**: an initial flank-length check
appeared to show zero flank added to every record. Root cause was in the
*check*, not the extraction: `bedtools getfasta` without `-nameOnly`
labels each output sequence with the coordinates of the BED interval it
was given (i.e. the already-flanked window), not the original core hit
coordinates — the first verification script assumed the header showed raw
coordinates and effectively double-subtracted the flank. Re-verified
properly by cross-referencing each alignment header's interval against the
real raw hit coordinates in `assigned.fasta` (not the alignment file's own
header) for 45 records across 9 files, 3 different subfamilies' up/down
extensions: all 45 match the confirmed boundary within 2bp. The extraction
was correct all along.

Published: 78 alignment files (29 top100 + 29 rand100 + 19 subfam — the
other 10 subfamilies either have <400 members, per the original
`report.html` stats table, or hit the getfasta bug above) replace the
previous fixed-50L/70R set entirely; `subfam/` rebuilt to match (19 `.al`
files, down from 25 since 6 more subfamilies' subfam variant is now
correctly skipped rather than silently broken). Full boundary numbers:
[`boundary_refine/boundary_refinement.tsv`](boundary_refine/boundary_refinement.tsv).

## Not yet done

- Re-verify/re-publish `eri`'s alignments with this same `(+,-)` fix (see
  above — real, not yet done).
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
