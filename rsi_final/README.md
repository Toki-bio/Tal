# rsi_final

The final *Rhinolophus sinicus* (rsi) run, ready for inspection and further use. Page with the summary and all open items: https://claude.ai/artifact/PfHq2MnKWxT4NcA76qJRkM

## What this is
One complete SINEderella run (assembly `GCF_057876585.1`) started from the **21 established consensuses** (r1-r10, MEG-RL/RS/T2/TR, seven composites: `consensuses.clean.fa`), followed by the flank scan (stages 1-7) and stage 8 on the same run. Server: therioserver `~/rhin/rsi_fresh2/run_20261002_112935` (flank scan in its `flankscan/`). Code: SINEderella `0f8f165` and later (mafft `--threadit 0` everywhere, `c34d052`; length-version test; consensus audit). Alignments are reproducible (same input, same alignment, any thread count); assignment counts vary by a few percent between runs because ssearch36 `-z 11` shuffles by design (final bank, length verdicts and the audit table are identical between two complete runs).

## Contents
| Path | What |
|---|---|
| `report.html` | full SINEderella report: assignment, divergence, plates (links to `alignments/`), element hierarchy, sequence similarity, length versions, consensus audit |
| `consensuses.clean.fa`, `consensuses.aln.fa` | the 21-consensus bank and its alignment (use as the bank of the next run, or with `--add` / `--exclude`) |
| `alignments/` | top 100, random 100 and subfamily plates for every unit and composite (60 files) |
| `stage8_plates/`, `stage8_tsd/` | single copies / copies inside composites / control (r9) plates and the TSD tables per family. Plates were built with `EXT3=1`: the proposed consensus row runs to the end of the conserved 3' tail (r5 +11 bp, r10 +21, r7 +8, others 0); element columns are unchanged, checked on all 18 plates |
| `flankscan/` | hierarchy, peaks, candidates, re-assignment and chain tables, accepted candidates (`cand/` alignments) |
| `length_variants/` | length-version test per pair (end modes, linkage, TSD); `summary.tsv` |
| `consensus_audit/` | each consensus rebuilt from its own copies (`summary.tsv`, `rebuilt.fa`) |
| `consensus_pairs/` | alignments of the pairs that matter for the length-version and subfamily questions |
| `r8_dimer_P42_P43_earlier_run.aln.fa` | combined alignment of the r8 + r8 dimer candidates from the earlier flank-scan run (in this run those copies are labelled with other composites) |
| `array_flag.tsv` | share of each family's copies in tandem arrays of regular spacing (SINEderella `tools/array_flag.py`): MEG-RS 92.8 % (21 arrays), every other family ~0 % |
| `assignment_stats.tsv`, `consensus_blocks.tsv` | per-subfamily assignment numbers; shared blocks between consensuses |

## Results in one place
* Length versions: r9 / r7, r9 / r8, r7 / r8 are separate SINEs (TWO_VERSIONS); r7 / r5 one SINE with a variable end; MEG-RS / MEG-RL and r4 / r3 not testable here (few copies).
* Audit: 11 MATCH, 5 SHORTER, 1 LONGER (r7 rebuilds 179 bp against 154: copies carry the r8 stretch), 1 DIVERGED (MEG-RS), 3 SKIPPED (MEG-RL, r2, r4).
* **MEG-RS is a tandem array in rsi** (1 592 of 1 715 copies in 21 arrays, median spacing 2 160 bp, flanks 90-92 % identical): the report says "Tandem array", the plates take independent copies first (top 100: 23 contigs, 32 rows marked `[array]`, still 2.8 kb wide). Not judgeable as a dispersed SINE from this genome.
* Flank scan on the bank with its composites: four new pairs of units pass the 70 % rule (r10 + P48, P48 + r8, C11 + P1, P26 + P34); nothing is added to the bank without your decision.

## Open (yours)
MEG-RS: keep as a tandem-array family or test on a bat genome where it is dispersed; r8 + r8 dimer; accept the four new pairs; peel subfamilies and rerun with `--add`; G-rich 3' end of P18; CpG-corrected divergence; r4 / r2 on another species.

## Known limits
* Plate step logs `border scan failed: No module named 'boundary'` (also in the original run): the border loop ran for 1 subfamily only; not yet investigated.
* `r3_58seqs_S` plate: one element column differs from the plain column majority (a tie, not an error).
* Earlier folders (`rsi/`, `rsi_v2` ... `rsi_v7`, `rsi_fresh`, `rsi_fresh2`) are superseded; `rsi_fresh2` is this run before the flank scan and plates.
