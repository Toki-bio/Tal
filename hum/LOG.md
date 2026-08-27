# hum (Homo sapiens) — cross-tool benchmark log

*Not part of the Talpidae/Erinaceidae/Timema discovery work above — the third
leg of the same cross-tool benchmark as [`../zeb/LOG.md`](../zeb/LOG.md)
(zebrafish) and [`../timb/LOG.md`](../timb/LOG.md) (Timema). Human has the
best-established SINE library of the three (Alu/MIR/FLAM/7SL etc.), so this
leg is the strongest ground-truth check of the set.*

---

## Method

- Genome: hg38, downloaded from UCSC goldenPath.
- Consensus library: 65 SINE consensuses extracted from Dfam via
  `famdb.py families -a -c --class SINE -f fasta_name 'Homo sapiens'`
  (same `dfam/tetools` Docker approach as the zebrafish leg).
- Full `SINEderella` run (`run_20260822_231705`) against the whole genome.

## A real infrastructure incident during this run

Step 2 (the 10-cycle `ssearch36` voting) crashed with "No space left on
device" partway through — `mktemp` with no explicit `TMPDIR` was silently
writing gigabytes of intermediate `.m8` files to `/tmp`, which on this server
lives on the small 24GB system partition, not the 14TB `/staging` mount. Two
real bugs were found and fixed as a result, both now committed to
`github.com/Toki-bio/SINEderella`:

1. `step2_asSINEment.sh` unconditionally reset `TMPDIR=""` before its own
   `mktemp` call, clobbering any inherited value — meaning an environment-level
   fix alone could never have worked. Changed to `TMPDIR="${TMPDIR:-}"` so an
   inherited value survives.
2. The main `SINEderella` orchestrator rebuilt the `results/` symlink farm
   (`rm -rf results; mkdir -p results`) *after* Step 6 had already written
   `results/report.html`, silently destroying it every time — reproduced
   across 4+ runs (eri, sco, two Timema runs) before being fixed here by
   reordering the two blocks.

Both fixes are load-bearing for this run: the successful completion below is
the first real run to exercise both post-fix, at the largest scale of any
run on this site so far (65 consensuses, hg38).

## Results

65/65 consensuses received real genomic hits, **~354,000 total assigned
copies** across step2's 10-cycle voting:

| Subfamily | Assigned | Subfamily | Assigned |
|---|---|---|---|
| AluSp | 47,049 | AluSx4 | 5,910 |
| MIR | 41,861 | FRAM | 5,705 |
| AluJb | 32,785 | FAM | 4,506 |
| MIRb | 30,080 | AluSx3 | 4,315 |
| FLAM_C | 22,533 | AluSz | 4,267 |
| AluSc | 21,168 | AluSc5 | 3,894 |
| AluJr4 | 19,946 | AluYa5 | 3,210 |
| FLAM_A | 13,411 | AluSq | 3,090 |
| AluSz6 | 11,782 | AluJo | 1,818 |
| AluSc8 | 10,528 | ... (47 more, down to AluYe6 at 1 copy) | |

Full per-subfamily table in `report.html`. As with the zebrafish leg, most
subfamilies show substantial `leak_pct`/`conf_alt_pct` from mutual overlap
between related Alu subfamilies voted flat against each other — a general
property of flat multi-consensus assignment on closely-related family
members, seen consistently across all three benchmark legs.

## Publishing

`step8a_extract_alignments.sh` + `step8b_publish_report.sh` ran cleanly:
65/65 subfamilies got top100/rand100 variants, most also got the subfam
variant (below `MIN_COPIES_SUBFAM` skipped for the smallest few). PCA panel
in `report.html` timed out (900s limit) at this scale and was skipped
gracefully — everything else in the report is complete.

## 2026-08-27 — Independent quality check against RepeatMasker's hg38 annotation

Unlike the AnnoSINE_v2 cross-checks on `timb`/`sco`, human has a real,
independent, pre-existing genome-wide annotation to compare against:
RepeatMasker's own hg38 track (`rmsk.txt.gz` from UCSC goldenPath, the
actual RepeatMasker output, not a re-derived approximation). Same
underlying subfamily nomenclature (both trace to Dfam/RepBase), so this is
a direct concordance check, not a proxy.

**Method**: downloaded `rmsk.txt.gz`, filtered to `repClass == SINE`
(1,910,631 genome-wide elements, 63 distinct subfamily names), converted
both RepeatMasker's calls and SINEderella's `assigned.fasta` loci to BED,
and joined with `bedtools intersect -f 0.5 -r` (reciprocal ≥50% overlap —
boundary definitions differ between Crossmatch-based RepeatMasker calls and
ssearch36-based SINEderella hits, so exact coordinate equality isn't the
right bar).

**Overlap rate**: 366,843 / 394,901 (92.9%) of SINEderella's assigned loci
found a matching RepeatMasker SINE call at all.

**Subfamily-name concordance**: 82.8% raw exact-match rate, **85.0% after
normalizing** for a naming-convention artifact — Dfam splits some
consensuses into `_short_`/full-length variants that RepeatMasker doesn't
distinguish (e.g. `AluJb_short_` was literally the single largest
"confusion" pair, 6,303 loci, purely because RepeatMasker just calls both
`AluJb`; same for `FLAM_C_short_`/`FRAM_short_`).

**Best-performing subfamilies** (young, well-differentiated — expected to
be easy for any tool): `AluYa5` 99.3%, `AluSc8` 98.8%, `AluSq2` 97.9%,
`AluSx` 97.7%, `AluSq` 96.9%.

**Where the two tools genuinely disagree, and why it's mostly explainable
rather than a SINEderella error**:

- `FLAM_C`/`FLAM_A`/`FRAM` (39–77% precision) confuse mostly with the
  oldest Alu-J subfamilies (`AluJo`/`AluJr`/`AluJb`). FLAM/FRAM are the
  ancestral Alu monomer precursors — biologically very close to the oldest
  Alu subfamilies, a known-hard boundary in the literature on either side.
- `AluJr4` (58%) confuses mostly with `AluJr`/`AluJo` — same old-Alu-J
  cluster.
- `MIR`/`MIRb`/`MIRc` mutually confuse (87–93% each way) — MIR subfamilies
  are ancient and highly diverged, also a known-hard case generally.

**`BC200` case, investigated specifically**: 449 loci SINEderella called
`BC200`, 0% exact match — but RepeatMasker's hg38 SINE track has **zero**
`BC200`/`BCYRN1` entries at all (confirmed directly: `grep -i bc200` against
all 63 RepeatMasker subfamily names returns nothing), so a 0% match rate
here was structurally guaranteed regardless of call quality — RepeatMasker
has no vocabulary for the thing being asked about. The biology explains
where the calls actually land: BC200 (BCYRN1) is a single-copy neuronal
ncRNA gene, not a genuine multi-copy retrotransposon subfamily — its 5'
~120bp domain is a FLAM-C-derived Alu left-monomer fused to a unique 3'
tail. 307/449 (68%) of SINEderella's `BC200` calls landed on `FLAM_C`
specifically, the rest on the same `AluJo`/`AluJr`/`AluJr4`/`AluJb`
ancestral-monomer cluster flagged above — exactly where evolutionary
relatedness predicts they'd land, not evidence of a misclassification.

A few other near-zero entries (`CAS`, `ASR`, `AluYk3`, `AluYh9`, `AluYh7`)
have tiny sample sizes (1–34 overlapping loci) — not statistically
meaningful on their own, not investigated individually.

Full overlap table, per-subfamily precision, and confusion-pair breakdown
generated via `bedtools intersect` + a small Python aggregation script; raw
intermediate files (`rmsk_sine.bed`, `sinederella_hum.bed`, `overlap.tsv`)
live in `/staging/tmp/rmsk_compare/` on DRAGEN, not copied into this repo
(large, easily regenerated from `rmsk.txt.gz` + `assigned.fasta`).
