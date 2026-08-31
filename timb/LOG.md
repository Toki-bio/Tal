# timb (Timema cristinae) — cross-tool benchmark log

*Not the curated Timema candidate-discovery work in [`../tim/`](../tim/) —
that project's consensuses come from this site's own manual `SubFam`
clustering. This is a separate benchmarking exercise: taking AnnoSINE_v2's
55 structurally-independent de novo candidates (found with no homology
library and no knowledge of `tim/`'s curated consensuses) and running each
one through SINEderella individually, purely to see how much genome each
pulls in and how much of the 55 is noise vs. real signal. See
[`../zeb/LOG.md`](../zeb/LOG.md) for the other benchmark leg using a
species with an established library instead.*

---

## Why this benchmark

Unlike zebrafish/human, Timema has no pre-existing SINE library, so this
benchmark can't compare against known-good consensuses. Instead it checks
something else: does an independent, purely-structural discovery tool
(AnnoSINE_v2 — TSD/poly-A/RNA-pol-III-promoter signal, no homology backstop
for this taxon) converge on the same subfamilies the manual curation in
`tim/` found, and how much noise sits in what it finds beyond that?

## Method

- **Server: DRAGEN** (`copilot@Biotech2024`, Tailscale `100.104.25.22`).
- **Run root**: `/staging/tmp/sinederella_benchmark/runs/timema/run_20260822_231705`
  — assignment output at `step2/step2_output/`, alignments at
  `results/alignments/`, report at `results/report.html`. Same genome file
  (`GCA_050494535.1.fa`) as `tim/`, different consensus set (AnnoSINE_v2's
  55 candidates, not the curated t1–t8).
- AnnoSINE_v2 mode 3 (hybrid) run directly on the Timema genome —
  55 candidates, all independently tagged tRNA-derived, zero hits from its
  built-in homology/HMM library (no related family exists in it for this
  taxon), so all 55 rest entirely on the structural module.
- Each of the 55 candidates run through `SINEderella` **individually as its
  own consensus** against the full genome (same tool, same assignment
  mechanism used everywhere else on this site) — not fed in together as one
  55-way vote.
- Cross-checked against `tim/`'s 6 manually-curated subfamilies
  (t1/t2/t345/t6/t7/t8, the set current at the time of this check) via
  `blastn` both directions (`-evalue 1e-10`, `-word_size 11`).

## Cross-check result

**All 6 curated subfamilies** have at least one significant AnnoSINE_v2 hit
(e-values 1e-12 to 1e-95, 82–89% identity). **40 of the 55 candidates (73%)**
independently converge onto one of the 6: t1 (14 matches), t2 (17), t345 (3),
t6 (6), t7 (1), t8 (2). Identity is consistently 76–89% — high enough to
clearly be the same element, low enough to read as two independent tools
agreeing, not one derived from the other.

## Noise breakdown on the 15 unmatched candidates

Matched vs. unmatched don't split cleanly by leak/sim_mean alone (matched
mean leak 73.2%, unmatched mean leak 66.3% — unmatched isn't uniformly
worse). Looking at individual candidates instead:

- **Genuinely clean, likely real missed signal**: `SINE_25` (15,940 copies,
  leak 0.18%), `SINE_47` (13,262 copies, leak 0.00%), `SINE_17` (260 copies,
  leak 0.00%) — as unambiguous as the best-behaved matched candidates,
  substantial copy counts, not noise.
- **Actually noisy**: `SINE_46`, `SINE_38`, `SINE_50`, `SINE_3`, `SINE_52`,
  `SINE_39`, `SINE_6`, `SINE_19`, `SINE_28`, `SINE_22`, `SINE_53`, `SINE_30` —
  small-to-modest copy counts with high leak/ambiguity, consistent with
  either fragments overlapping existing families too closely to cleanly
  match, or structural false positives from AnnoSINE_v2's purely-structural
  (no homology-backstop) calling on this taxon.

## A real bug found and fixed here: reverse-complement orientation

Several alignments (e.g. `SINE_18` top100) showed roughly half the sequences
genuinely reverse-complemented relative to the consensus — confirmed by
k-mer-vs-consensus orientation scoring (48/100 sequences matched the
consensus's reverse complement better than the consensus itself), despite
every header claiming `(+)` strand. This is isolated to this benchmark run —
the curated `tim/` publish checked clean (100/100 forward) — strand tracking
broke somewhere upstream specifically for this externally-sourced candidate
bank, not a general pipeline defect. Fixed by re-aligning all 144 alignment
files with `mafft --adjustdirection`, which auto-detects and flips
misoriented sequences per profile; re-verified via the same k-mer check
(`SINE_18`: 48 reverse → 0 reverse after the fix). `step8a`'s default
alignment step does not use `--adjustdirection` — worth adding when a
consensus bank's strand-tracking provenance isn't guaranteed clean (e.g. any
externally-sourced candidate set like this one), flagged as a follow-up.

## Publishing

`report.html` is the standard SINEderella-generated report (same generator
as every other species page). All 55 candidates got top100/rand100
alignments; 44/55 also got the subfam variant (11 below the
`MIN_COPIES_SUBFAM` threshold).
