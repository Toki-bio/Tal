# zeb (Danio rerio) — cross-tool benchmark log

*Not part of the Talpidae/Erinaceidae/Timema discovery work above — this is a
separate benchmarking exercise, comparing SINEderella against AnnoSINE_v2 on a
species with an established SINE consensus library, so the "ground truth" is
known going in. See [`../timb/LOG.md`](../timb/LOG.md) and the human run for
the other two legs of the same three-species comparison.*

---

## Why zebrafish

Zebrafish (*Danio rerio*, danRer11/GRCz11) has a well-curated SINE library in
Dfam (SINE2/HE/SINE3/DANA/7SLRNA families etc.), unlike Timema which has none.
Running SINEderella against a genuinely independent, pre-existing consensus
set — rather than consensuses derived from the same discovery run being
tested — makes this a real benchmark, not a circular check.

## Method

- **Server: DRAGEN** (`copilot@Biotech2024`, Tailscale `100.104.25.22`).
- **Run root**: `/staging/tmp/sinederella_benchmark/runs/zebrafish/run_20260822_231705`
  — assignment output at `step2/step2_output/`, alignments at
  `results/alignments/`, report at `results/report.html`.
- Genome: danRer11 (GRCz11), downloaded from UCSC goldenPath.
- Consensus library: 15 SINE consensuses extracted from Dfam via
  `famdb.py families -a -c --class SINE -f fasta_name 'Danio rerio'`
  (the `dfam/tetools` Docker image; the actual `.h5` family database lives at
  `Libraries/famdb/` inside the image, not the top-level path the CLI's own
  `-i` flag description implies — a real gotcha hit getting this working).
- Full `SINEderella` run (`run_20260822_231705`) against the whole genome.

## Results

15/15 consensuses received real genomic hits, **155,249 total assigned
copies**:

| Subfamily | Total assigned | Firm % | Leak % | Sim mean |
|---|---|---|---|---|
| HE1_DR1 | 64,343 | 37.43 | 31.79 | 0.676 |
| SINE2-3_DR | 33,298 | 20.58 | 98.97 | 0.778 |
| HE2_DR | 13,658 | 8.53 | 98.58 | 0.936 |
| SINE2-4_DR | 12,152 | 4.55 | 64.61 | 0.659 |
| SINE2-5_DR | 6,570 | 2.64 | 32.80 | 0.661 |
| SINE3-1a | 4,049 | 2.55 | 0.03 | 0.726 |
| SINE2-1_DR | 4,461 | 2.41 | 0.13 | 0.606 |
| SINE3-1 | 3,490 | 2.21 | 0.03 | 0.775 |
| SINE2-1B_DR | 1,841 | 0.92 | 0.21 | 0.604 |
| SINE_DR1 | 1,874 | 0.60 | 0.00 | 0.519 |
| SINE_TE | 1,552 | 0.28 | 86.84 | 0.536 |
| 7SLRNA | 212 | 0.14 | 0.00 | 0.552 |
| SINE2-2_DR | 236 | 0.12 | 0.00 | 0.406 |
| LmeSINE1c | 446 | 0.10 | 49.07 | 0.410 |
| DANA | 7,067 | 0.03 | 93.02 | 0.396 |

Most subfamilies show high `leak_pct` (many >60-95%) — a lot of mutual
ambiguity between these 15 Dfam consensuses when voted on flat and unweighted
against each other, similar to what's seen benchmarking AnnoSINE_v2's 55
Timema candidates (`../timb/LOG.md`) — this looks like a general property of
flat multi-consensus voting with related/overlapping families, not specific
to either species or tool.

## Publishing

`step8a_extract_alignments.sh` + `step8b_publish_report.sh` ran cleanly.
Initially 13/15 subfamilies got all three alignment variants
(top100/rand100/subfam); the other 2 (`7SLRNA`, `DANA`, `LmeSINE1c`,
`SINE2-2_DR` — small counts) were correctly skipped for `subfam` below the
`MIN_COPIES_SUBFAM` threshold.

**A real bug found here**: 3 subfamilies with substantial copy counts
(`SINE2-3_DR` 32,619, `SINE2-4_DR` 7,267, `SINE3-1a` 4,047) were *also*
missing their `subfam` variant, and not because of the threshold —
`bedtools getfasta` was failing outright with "malformed BED entry, Start
Coordinate detected that is < 0". Root cause: the BED-construction step
(`$4-1` for 0-based conversion) can produce a negative start when a locus
sits right at the beginning of a contig, and `bedtools getfasta` aborts its
*entire* batch on a single malformed line. `top100`/`rand100` (100 loci
sampled each) rarely hit this by chance; `subfam` samples up to 10,000 loci,
so the odds of including one contig-edge locus are much higher — explaining
why only the `subfam` variant was silently and completely lost for these
three, while `top100`/`rand100` for the same subfamilies worked fine.
Checked `tim`/`timb`/`hum`'s manifests for the same pattern (`has_subfam=0`
despite ≥400 members) — none affected, so this appears to be specific to
where zebrafish's SINE loci happen to fall relative to contig boundaries in
this particular assembly. Fixed by clamping the start coordinate to 0 in
`step8a_extract_alignments.sh` (both occurrences of the BED-construction
awk); re-ran the `subfam` extraction for the 3 affected subfamilies with the
fix, all three now have the full three-variant set. Fix committed to
`github.com/Toki-bio/SINEderella`.
