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
