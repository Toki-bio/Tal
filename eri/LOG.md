# eri (Erinaceus europaeus) — analysis log

*Scoped to this species only — see the main [WORKFLOW.md](../WORKFLOW.md) for the
shared Tal SINE pipeline this reuses. Genome: GCF_950295315.1 (mEriEur2.1).*

---

## 2026-08-20 — Initial run: hit a real SubFam bug

First `SINEderella` invocation on this genome
(`GCF_950295315.1_mEriEur2.1_genomic.fna`, 2 known consensuses ERI-1/ERI-2 from
Borodulina 2001, run directory `run_20260820_204827` on KIT) completed genome
sanitization, ran Step 1 search & extract successfully (1,548,346 candidate loci
extracted, 251,355 hits for ERI-1 / 1,312,311 hits for ERI-2), then **failed** at
Step 1's SubFam clustering sub-step:

```
/data/V/toki/bin/SINEderella/SubFam: 7: /data/V/toki/bin/SINEderella/SubFam: [[: not found
/data/V/toki/bin/SINEderella/SubFam: 11: /data/V/toki/bin/SINEderella/SubFam: Syntax error: "(" unexpected
```

**Root cause**: `SubFam`'s shebang was `#!/bin/sh`, but the script body used
bash-only syntax (`[[ ]]` test brackets, array assignment). On KIT, `/bin/sh` is
dash, which doesn't support either — 100% reproducible on every invocation.
Confirmed identical bug in the canonical GitHub copy
(`Toki-bio/SINEderella`), not local drift. Also confirmed the same failure had
already silently occurred in an earlier, unrelated run (the dog/*Canis
familiaris* validation run) without being noticed, since that run's Step 2
onward was resumed manually, bypassing the orchestrator's fatal-exit check on
Step 1.

**Fixed**: shebang changed to `#!/bin/bash` (single-line diff). Deployed to the
live shared tool path (`/data/V/toki/bin/SINEderella/SubFam`, PATH-resolved by
every invocation, not a per-run frozen copy) and committed/pushed to
`Toki-bio/SINEderella` (`main`, commit `8feb497`/merge `8c65436`), alongside two
other fixes made the same session (LEAK-detection dead-code fix in
`step2_asSINEment.sh`, python3-version auto-detect fix in `step4_plots.sh`).

A second attempt to restart the full `SINEderella` pipeline was interrupted by
the user mid-genome-sanitization (`run_20260820_210902`) — that run's
`genome.clean.fa` was left at 0 bytes (confirmed directly), so it was abandoned
rather than resumed (SINEderella's `--resume` mode only checks file
existence/non-emptiness, not completeness, so resuming a truncated genome file
would have silently run the full pipeline on corrupt input).

## 2026-08-20 — Reused Step 1's output, ran Step 2 manually

Rather than redo the expensive Step 1 search from scratch, reused the *original*
run's already-complete `genome.clean_step1/extracted.fasta` (1,548,346 sequences,
produced before SubFam killed the run) directly as Step 2's input — same pattern
used for the dog genome validation run. Deployed the fixed (LEAK-tracking)
`step2_asSINEment.sh` into `run_20260820_204827/`, launched:

```bash
bash step2_asSINEment.sh consensuses.clean.fa genome.clean_step1/extracted.fasta 48 step2/step2_output
```

**Result**: 1,548,011/1,548,346 unanimous (10/10 votes), 1,276,402 passing all
criteria (82.4%) — ERI-1: 236,299 sequences, ERI-2: 1,311,712 sequences.

## 2026-08-20 — Ran the now-fixed SubFam directly against the original sample

`genome.clean_step1/subfam_input/input.fasta` (the 30,000-sequence sample SubFam
was originally given before crashing) still existed from the first attempt.
Ran the fixed `SubFam` directly against it (not through the full orchestrator,
since Step 1 as a whole was bypassed):

```bash
export PATH=/data/V/toki/bin/SINEderella:$PATH
cd run_20260820_204827/genome.clean_step1/subfam_input
SubFam input.fasta 50
```

Completed successfully in 4m9s ("Reordering bank and making consensus
sequences" → "Aligning consensus sequences" → "Final alignment" → "Job
completed"), producing `input.clw` — an alignment of ~600 per-chunk consensus
sequences (each chunk MAFFT-aligned independently, then all chunk-consensuses
aligned together in a final cross-chunk step), not a raw alignment of the
original 30,000 sequences and not true CLUSTAL block format despite the `.clw`
extension (plain gapped FASTA). Downloaded locally.

## 2026-08-20/21 — Manual clustering into candidate subfamilies

User manually inspected `input.clw` and clustered the per-chunk consensus
sequences into two candidate subfamily groups:
- **ERI-1-like**: e1-1 (36 seqs, closest match to the literature ERI-1
  consensus), e1-2 (19 seqs), e1-3 (32 seqs), e1-4 (3 seqs)
- **ERI-2-like**: e2-1 (225 seqs), e2-2 (4 seqs), e2-3 (231 seqs), e2-4 (44 seqs)

Result aligned against the original SINEbase/literature ERI-1 and ERI-2
consensus sequences (Borodulina 2001) — see
[`alignments/eri_candidate_consensuses_vs_sinebase.aln.fa`](alignments/eri_candidate_consensuses_vs_sinebase.aln.fa).

**Next step (user, not yet done)**: run `SINEderella` against this new,
expanded candidate consensus set (8 sequences instead of 2) to get a real
subfamily classification, then extract top100/rand100/subfam alignments per
subfamily (see [`scripts/`](scripts/)) and build a full `report.html` matching
the other 6 Tal species.

## 2026-08-21 — Second SINEderella run against the 8-consensus candidate set

Ran the full pipeline (`run_20260820_221537` on KIT) using
`eri_candidates8.fasta` → `consensuses.clean.fa` (8 sequences: e1-1..e1-4,
e2-1..e2-4) as the reference bank. Step 1 (search+extract) and Step 2
(assignment) both completed. All 8 candidates received real genomic hits —
none came back empty:

| Subfamily | Assigned copies |
|---|---|
| e1-1 | 123,749 |
| e1-2 | 146,021 |
| e1-3 *(header has no space: `e1-3consensus_32seqs`, inherited from the source consensus file, not fixed here)* | 83,033 |
| e1-4 | 10,859 |
| e2-1 | 631,692 |
| e2-2 | 11,578 |
| e2-3 | 647,525 |
| e2-4 | 81,431 |

## 2026-08-21 — Extraction script: found and adopted the real Tal precedent

The first version of `scripts/extract_top100_rand100_subfam.sh` (staged
2026-08-20, described above) was written without an exact precedent —
its own header said so explicitly. Two problems surfaced when actually
running it against real data:

1. **Wrong sampling method.** It used a custom seeded Python
   `random.sample()` for `rand100`. The project has its own dedicated tool
   for this, `/data/V/toki/bin/sample N file.fa` (unseeded, `shuf`-based),
   which is the one that should be used for reproducing the same convention
   as the rest of the site.
2. **Missing genomic flanks.** The first version aligned the raw matched
   region only. The real Tal species pages include 50bp-left/70bp-right
   genomic flanks around each copy (for TSD inspection) — confirmed by two
   independent pieces of evidence: the Tal repo's own commit history
   (`"Add flanked alignments (50+70bp)"`, `"...flanked versions (50bp left,
   70bp right)"`) and direct user confirmation.

Searched KIT directly for the actual script that produced the existing
saq/ccr/etc. alignments (git history in the Tal repo itself only ever
contains the *output* `.aln.fa` files, never a generation script — so this
had to be found on the server, not in either git repo). Found
`/data/V/toki/bin/SINEderella/extract_alignments.sh`, confirmed identical
(md5-matched, modulo text encoding) to `/data/W/toki/
extract_alignments_sq.sh` — the actual script invoked for a prior Tal
run (log at `/data/W/toki/extract_alignments_sq.log`). It implements exactly
this: re-extract flanked copies from the genome via `bedtools slop`+
`getfasta` using coordinates parsed from the SINEderella header format
(`>ctg:start-end(strand)|Subfam|Bits`), then a post-alignment cleanup step
that lowercases and re-justifies the flank bases against the consensus body
boundary. Its own header comment says "30bp left", which conflicts with the
50bp confirmed by the commit history/user — used 50, not 30.

`scripts/extract_top100_rand100_subfam.sh` was rewritten to fold in this
exact real logic (functions renamed/kept close to the original for
traceability) instead of the original from-scratch reimplementation. Key
parts, embedded here since they're the actual mechanism behind every
`top100`/`rand100` alignment on this page:

```bash
# Parse FASTA headers to BED. Header format: >ctg:start-end(strand)|Subfam|Bits
fasta_headers_to_bed(){
  local infile="$1" outfile="$2"
  grep "^>" "$infile" | sed 's/^>//' | awk '{
    id = $1
    n = split(id, parts, "|")
    loc = parts[1]
    strand = "+"
    if (sub(/\([+-]\)$/, "", loc)) {
      i = index(parts[1], "(")
      strand = substr(parts[1], i+1, 1)
    }
    last_colon = 0; rest = loc
    while ((p = index(rest, ":")) > 0) { last_colon += p; rest = substr(rest, p+1) }
    if (last_colon == 0) next
    ctg = substr(loc, 1, last_colon-1)
    coords = substr(loc, last_colon+1)
    split(coords, c, "-")
    if (c[1]+0 >= 0 && c[2]+0 > 0) {
      printf "%s\t%s\t%s\t%s\t0\t%s\n", ctg, c[1], c[2], id, strand
    }
  }' > "$outfile"
}

# Re-extract with FLANK_L=50 / FLANK_R=70 flanks from the genome
extract_flanked(){
  local infile="$1" outfile="$2" tmpd="$3"
  fasta_headers_to_bed "$infile" "$tmpd/flank.bed"
  [[ -s "$tmpd/flank.bed" ]] || { cp "$infile" "$outfile"; return; }
  awk -v OFS='\t' '{print $1,$2}' "${GENOME}.fai" > "$tmpd/genome.sizes"
  bedtools slop -s -l "$FLANK_L" -r "$FLANK_R" -g "$tmpd/genome.sizes" \
    -i "$tmpd/flank.bed" > "$tmpd/flanked.bed" 2>/dev/null || true
  [[ -s "$tmpd/flanked.bed" ]] && bedtools getfasta -s -nameOnly -fi "$GENOME" \
    -bed "$tmpd/flanked.bed" > "$outfile" 2>/dev/null || true
  [[ -s "$outfile" ]] || cp "$infile" "$outfile"
}
```

Full script (including the `postprocess_flanks` alignment-cleanup step, also
carried over from the real precedent): [`scripts/extract_top100_rand100_subfam.sh`](scripts/extract_top100_rand100_subfam.sh).

Re-ran against `run_20260820_221537`; all 8 subfamilies completed. Built and
pushed `report.html` (via `step6_report.py --tal-species-code eri
--no-sineplot`, using this run's own `manifest.txt` — all stdlib, no extra
deps needed since SINEplot's own PCA panel is separate from the internal
mutation-landscape PCA `step6_report.py` builds itself, which ran fine) and
updated `index.html`'s eri card to "full report".

## 2026-08-21 — Two real bugs found in the published alignments after the fact

User caught by eye that the MSA-viewer view of `eri_e1-1_rand100.aln.fa`
showed sequences ending right at the SINE body with no visible flank — I
initially answered wrong, claiming the flank data was present but just not
visually distinguished by MSA-viewer's case-insensitive coloring (partially
true as a separate finding, see below, but not the actual explanation here).
Checked properly and found two real bugs:

**Bug 1 — no real flanks were extracted.** Direct measurement (ungapped
sequence length vs. core-region length parsed from the header coordinates)
showed `ungapped_length == core_length` for every sequence checked — zero
flank bases actually added, despite the script's own log printing
`"OK top100 (100 copies, 50L+70R flanked, + consensus)"` for every
subfamily. Root cause of the false confidence: that log line is printed
*unconditionally* after calling `extract_flanked`, regardless of whether its
internal `bedtools`-failure fallback (silently `cp`s the unflanked input)
triggered. `postprocess_flanks` then lowercased ordinary MAFFT-alignment
edge/indel overhang and made it *look* like real flank data (case-wise),
which is what fooled the first check. The `extract_flanked` mechanism itself
tested correct in every isolated and fresh-debug-run reproduction — the
original failure mode was never conclusively root-caused, so the fix was a
clean re-run with real per-record verification (see below) rather than a
trust-the-log check, plus removing reliance on the log message as any kind
of confirmation.

**Bug 2 — contig names leaked `@U@` instead of `_`.** E.g.
`NC@U@080165.1` instead of `NC_080165.1`. Found via
`grep -n '@U@' /data/V/toki/bin/SINEderella/SINEderella` →
`gsub(/_/,"@U@",hdr)` (line 214) — the orchestrator sanitizes underscores in
contig names for its own internal delimiter safety. `sear`/`sear_multi`
already restore `@U@`→`_` on their own outputs (confirmed:
`sear:588-607`, `sear_multi:378-382`), but this project's extraction
scripts never did, so the substitution leaked all the way to the published
alignment headers.

**Fix**: both `scripts/extract_top100_rand100_subfam.sh` and
`scripts/extract_tribes.sh` now (1) restore `_` from `@U@` on every output
file right before writing it, and (2) get their flanking verified with
direct measurement, not log trust — for the corrected re-run, 20 real
records sampled across 4 different output files were checked
programmatically: `ungapped_length == (header_core_length + 120)` for every
one, and a `grep -c '@U@'` sweep across all output files returned zero
matches, before any file was downloaded or committed.

**Separate, real finding**: MSA-viewer itself renders lowercase and
uppercase bases identically (`script.js:6301-6302`,
`getResidueAnnotationClasses()` calls `.toUpperCase()` before deriving any
CSS class) — so even with genuinely correct flank data, there is currently
no visual way to distinguish flank (lowercase) from body (uppercase) bases
in the viewer. Not yet fixed; would need a `flank` CSS class keyed off
original case before the uppercase-normalization step.

## 2026-08-21 — Tandem-repeat contamination: quantified genome-wide, checked per-Tribe

Follow-up on the periodicity finding above. Ran genome-wide Tandem Repeats
Finder (`/data/V/toki/bin/trf`, parallelized 32-way by contig via
[`run_trf_genome.sh`](scripts/run_trf_genome.sh) — not committed to
`scripts/` yet, see below) against `genome.clean.fa`, then intersected the
result against both the full SINE hit set and each Tribe individually.

**Two real bugs hit and fixed along the way** (documented here since they're
easy to repeat):
1. `trf` (the installed v4.09 build) segfaults above `MaxPeriod≈2000` — every
   value tried from 2200 up crashed it consistently. Empirically verified
   2000 is still sufficient: a direct test on tribe10's known periodic
   region found real tandem content at MaxPeriod=2000 (small microsatellite
   periods of 1/3/4/10bp, recurring in an alternating ~1,160bp/~2,260bp
   pattern that sums to the same ~3,425bp macro-period found earlier) — the
   true tandem unit is tiny, the larger inter-hit spacing is just the
   SINE-query-motif's recurrence interval within the array.
2. First full run "completed" in 1 second with 0 results — `$OUTDIR` was
   passed as a relative path, and each parallel job does `cd "$OUTDIR/dat"`
   before referencing a path built from `$OUTDIR`, breaking resolution
   after the cd. Fixed by resolving `$OUTDIR`/`$GENOME` to absolute paths
   upfront. Second bug, same symptom class: `trf2bed` writes its own
   `<input>.bed` file rather than printing to stdout (looks like it should
   from the command name, isn't) — the script was capturing stdout, getting
   nothing, silently succeeding with an empty merged BED. Fixed by `cat`-ing
   the real per-chunk `.bed` files it produces.

**Genome-wide result**: 1,105,283 merged tandem-repeat intervals, ~103.7Mb
total (~4% of the 2.6Gb genome).

**SINE hit exclusion estimate**: intersecting against all 1,735,888 assigned
SINE hits (`step2_output/assigned.fasta`), **207,492 (11.95%) fall inside a
tandem-repeat interval** and should be excluded from the true SINE copy
count. Per subfamily:

| Subfamily | Total | In tandem | % |
|---|---|---|---|
| e1-1 | 123,749 | 13,244 | 10.7% |
| e1-2 | 146,021 | 29,917 | 20.5% |
| e1-3 *(header: e1-3consensus_32seqs)* | 83,033 | 6,521 | 7.9% |
| e1-4 | 10,859 | 85 | 0.8% |
| e2-1 | 631,692 | 74,333 | 11.8% |
| e2-2 | 11,578 | 141 | 1.2% |
| e2-3 | 647,525 | 75,876 | 11.7% |
| e2-4 | 81,431 | 7,375 | 9.1% |

**Per-Tribe check — finding "proper" tribes with unique flanks**: this is
the direct, authoritative check (vs. the periodicity-in-hit-spacing proxy
used in the previous entry, which turned out to be an unreliable predictor —
see below).

| Tribe | Tandem overlap (of 500 sampled members) |
|---|---|
| tribe02 | **89.0%** — contaminated |
| tribe07 | **89.6%** — contaminated |
| tribe09 | **88.6%** — contaminated |
| tribe01 | 0.2% — clean |
| tribe03 | 0.2% — clean |
| tribe04 | 0.0% — clean |
| tribe05 | 0.0% — clean |
| tribe06 | 0.2% — clean |
| tribe08 | 0.0% — clean |
| tribe10 | 0.0% — clean |

**7 of 10 tribes (01, 03, 04, 05, 06, 08, 10) are genuinely clean, dispersed
SINE families** — not tandem-array artifacts.

**The periodicity proxy from the previous entry was a poor predictor**: it
flagged tribe06/08/09/10 as suspicious (strong periodic SINE-hit spacing,
59-94% regularity) and called tribe02 the cleanest (weakest periodicity,
19%). Direct TRF evidence agrees only on tribe09; it's wrong for tribe02
(actually 89% contaminated despite irregular spacing) and wrong for
tribe06/08/10 (actually clean despite regular spacing). Lesson: periodic
hit-*spacing* and actual tandem-repeat *sequence content* overlap are not
the same signal — a tribe can sit mostly inside satellite DNA with
irregular hit spacing, or have regular hit spacing while sitting almost
entirely outside any TRF-detected tandem region. Direct intersection
against real tandem annotation (this entry) is the reliable check; the
spacing heuristic (previous entry) was only useful for raising the
question, not answering it.

Scripts: [`scripts/run_trf_genome.sh`](scripts/run_trf_genome.sh),
[`scripts/per_subfamily_tandem_overlap.sh`](scripts/per_subfamily_tandem_overlap.sh),
[`scripts/per_tribe_tandem_overlap.py`](scripts/per_tribe_tandem_overlap.py).
Data: [`tandem_analysis/summary.tsv`](tandem_analysis/summary.tsv) (the two
tables above). The full interval/hit BED files (`tandem_repeats.merged.bed`,
1.1M intervals; `sine_hits.in_tandem.bed`, the 207,492 excluded hits) are
*not* committed — repo convention excludes all `*.gz`/raw per-hit BED data
(see `.gitignore`, matches how every other species' `assigned.fasta`/
`all_hits.labeled.bed` etc. are handled). They remain on KIT in
`run_20260820_221537/` (`trf_out/tandem_repeats.merged.bed`,
`sine_hits.in_tandem.bed`) for anyone who wants to audit specific loci.

**Not yet done** (at the time of that entry): the 7 clean tribes hadn't been
re-published as a filtered/labeled subset — see next entry, superseded by
checking further down the size ranking instead.

## 2026-08-21 — Went past the top-10 clusters, found 11 genuinely clean ones

Re-ran the vsearch clustering (`--cluster_fast --id 0.99 --strand plus` on
`genome.clean_step1/extracted.fasta`) to regenerate `clusters.uc` (deleted
by the original run's cleanup — see the intermediates-preservation fix,
same commit that removed the auto-`rm -rf $WORK` from both extraction
scripts). Reproduced the original top-10 tribes exactly (same cluster IDs,
same sizes) — confirms the clustering is deterministic, not an artifact of
run-to-run variation.

Checked 50 more clusters, ranks 11&ndash;60 by size (133&ndash;373 members,
below the size range the original 10 "Tribes" were drawn from), for the
same two criteria as the tandem-repeat entry above: genuine per-member
flank uniqueness (singleton left+right flank, not shared with any other
cluster member) and TRF tandem-repeat overlap.

**11 of 50 pass both convincingly**:

| Cluster | Size | Singleton-flank | Tandem overlap |
|---|---|---|---|
| 449354 | 373 | 59.8% | 1.9% |
| 11535 | 229 | 79.0% | 12.7% |
| 159155 | 222 | 81.5% | 26.6% |
| 24069 | 199 | 72.4% | 13.1% |
| 10605 | 146 | 81.5% | 14.4% |
| 19132 | 159 | 74.2% | 1.3% |
| 189222 | 149 | 75.8% | 10.7% |
| 11960 | 149 | 70.5% | 20.1% |
| 447956 | 155 | 63.9% | 20.6% |
| 184652 | 133 | 66.2% | 4.5% |
| 158992 | 133 | 64.7% | 3.8% |

Manually spot-checked the two highest-singleton clusters (`11535`,
`449354`) by eye: left flanks are genuinely distinct strings across
different scaffolds at wildly different coordinates (e.g. `449354`:
`NC_080162.1` at 187M, 190M, 194M, 195M — real dispersed positions, not a
repeated motif). Right flanks in `11535` show a family resemblance
(AT-rich, some shared short substrings) but with real accumulated
divergence between copies, not byte-identical duplication — consistent
with genuine evolutionary drift from a common ancestral insertion, which
is exactly what real dispersed SINE copies are expected to look like.

The other 39/50 checked show the same contamination signature as the
original top-10 (0&ndash;1.8% singleton-flank, mostly exactly 0%) — so this
isn't "everything below the top 10 is clean," specifically these 11 stand
out.

Published as a new page:
[`eri/unique_flank_tribes.html`](unique_flank_tribes.html) (linked from
`index.html` and from `report.html`'s Tribes section) — real MAFFT
alignments (50L+70R flanks, same convention as the rest of the site) built
from each cluster's full membership, plus
[`eri/unique_flank_tribes/summary.tsv`](unique_flank_tribes/summary.tsv)
with the full numbers for all 50 checked (not just the 11 that passed).

**Not yet done** (at the time of that entry): hadn't checked clusters below
rank 60 — see next entry, and the "11 verified" claim below is retracted.

## 2026-08-21 — Correction: the "11 verified" claim doesn't survive a real statistical test

User caught the real problem with the previous entry's methodology: the
"singleton flank" test only checked for **exact-duplicate** flanks. That's
too weak — it doesn't catch flanks that are similar-but-not-identical,
which is exactly the signature ordinary mutation accumulation between
related copies produces. Direct check on `11535` (one of the "best" 11):
mean pairwise difference across its right flanks was only 33.5% (66.5%
identity) with a clearly recognizable shared motif across samples
(`TAACTACAACAATAAA...CAACAAGGGCAACAAAAGGGAATAAATAAAT...`) — obviously
related sequence, not independent DNA, that the exact-match test completely
missed since no two copies were byte-identical.

**Established the real background first**: sampled 2,000 pairs from 300
random unrelated 70bp genomic windows (all contigs) → **25.0% mean pairwise
identity** — matches the theoretical value for independent 4-letter DNA
exactly, confirming no meaningful compositional bias needs correcting for
in this genome.

**Re-scanned properly**: required both left AND right flank mean pairwise
identity (150 sampled member pairs each) to be ≤35% (background + margin)
to count as genuinely unique — not "no exact duplicate," but "statistically
indistinguishable from unrelated DNA." Scanned 48 clusters, ranks 11–500 by
size (all 11 originally-flagged clusters included).

**Result: 0 of 48 pass on both sides.** But a real, informative pattern in
the failures: left (upstream) flanks are often close to background (31–33%
in the previously-published 11) while right (downstream) flanks are
consistently far above it (58–100%, no exception, even in the "best"
clusters). Leading hypothesis, not yet confirmed: SINE 3' ends are
notoriously fuzzy (poly-A tails, AT-rich trailing sequence), so the
original de novo hit-boundary calling likely under-calls the true 3' end —
meaning the "downstream flank" window is still real SINE-tail sequence for
many copies, not independent genomic context. Not tested: whether shifting
the right-flank window further downstream (past the likely fuzzy boundary
zone) recovers background-level identity.

**Standing conclusion, corrected**: no cluster checked so far — the
original top-10 Tribes, nor these 48 further down the ranking — has been
shown to have genuinely independent flanking DNA on both sides. The "11
verified unique-flank clusters" claim from the previous entry is retracted.
[`eri/unique_flank_tribes.html`](unique_flank_tribes.html) updated in place
with a correction notice at the top, the original (now-labeled-retracted)
writeup below it, and the full 48-cluster strict-scan table with real
numbers
([`unique_flank_tribes/strict_scan_summary.tsv`](unique_flank_tribes/strict_scan_summary.tsv)).
The 11 alignment files are kept for reference/audit, not deleted, but
should not be read as verified real dispersed SINE families.

**Not yet done**: the fuzzy-3'-boundary hypothesis is untested. Also
haven't checked clusters below rank 500 (3,484 total clusters with ≥5
members exist) under the corrected criterion.

## 2026-08-21 — e2-3 boundary extension: the fuzzy-3'-boundary hypothesis, tested

User noticed e2-3 copies sometimes appear to continue past the called
consensus border (right side, occasionally left), and asked for a direct
test on e2-3's top100/rand100 sets rather than the Tribes clusters:
stepwise `bedtools slop` extension (50bp increments, cap 1000bp) on each
side, checking at each step whether the newly-added window's pairwise
identity across the 100-member sample has dropped to this genome's
established background (25.0%, +10 margin = 35% threshold — same numbers as
the Unique-Flank Tribes work above).

**Real bug hit and fixed while building this** (worth recording — same bug
class flagged repeatedly in other scripts this project): the script
restored `@U@` → `_` in contig names immediately after parsing
`assigned.fasta`, before using those names to query the genome via
`samtools faidx`. But `genome.clean.fa` itself keeps `@U@`-sanitized contig
names (confirmed directly: `genome.clean.fa.fai` lists `NC@U@080162.1`, not
`NC_080162.1`) — restoring underscores that early broke every coordinate
lookup, and the script silently produced empty/truncated output with no
error under `set -euo pipefail`, plus separately looked like a hang because
the SSH tunnel to KIT was dropping around the same time (both problems were
real, independently, and had to be untangled one at a time). Fixed by
deferring `@U@` restoration to only the final alignment output files.

**Result**:

| Set | Side | Extension needed | Identity at that point |
|---|---|---|---|
| top100 | upstream (5') | 0bp | 17.1% |
| top100 | downstream (3') | 50bp | 33.1% |
| rand100 | upstream (5') | 0bp | 14.8% |
| rand100 | downstream (3') | 50bp | 28.2% |

The upstream (left) side is already at background identity at the
originally-called boundary — no extension needed, the left boundary call
looks correct. The downstream (right) side needed a 50bp extension before
identity dropped to background in both independent samples (top100 and
rand100 agree) — **this confirms the fuzzy-3'-boundary hypothesis** raised
in the previous entry: the original de novo 3' boundary call under-calls
the true SINE end by roughly 50bp for e2-3.

Published as a new section on
[`report.html`](report.html#e2-3-extended) right after the main alignments
table:
[`alignments/eri_e2-3_top100_extended.aln.fa`](alignments/eri_e2-3_top100_extended.aln.fa)
and
[`alignments/eri_e2-3_rand100_extended.aln.fa`](alignments/eri_e2-3_rand100_extended.aln.fa)
(50bp upstream + 120bp downstream flanks — base 50/70 convention plus the
confirmed 50bp extension on the downstream side), raw boundary-scan numbers
in [`e2-3_extend/boundary_results.tsv`](e2-3_extend/boundary_results.tsv).

**Not yet done**: only e2-3 has been tested this way; e1-1..e1-4, e2-1,
e2-2, e2-4 have not been checked for the same fuzzy-boundary effect. This
same stepwise-extension method is also being generalized into a proper
SINEderella pipeline step (`step7_boundary_refine.sh`, applied per-locus
rather than per-sample-set) — not yet run against real eri/scorpion data as
of this entry.
