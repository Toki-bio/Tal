# rsi_main_20261018

*Rhinolophus sinicus* (rsi), complete SINEderella run with **flank verdicts**: same genome (`GCF_057876585.1`), same
21-consensus bank and settings as `rsi_main_20261016`; only the code changed.

- Code: SINEderella `main` 6411da2 (2026-10-18). New against 296fa33: the flank expansion's end state decides
  (`tools/flank_verdict.py`), every copy gets a verdict, and the plates show SINE loci only (`tools/plate_filter.py`).
  Aim of the workflow: SINEderella `docs/AIM.md` (interspersed repeats = copies with unique flanks).
- Server: therioserver `/home/toki/rhin/rsi_main_20261018/run_20261017_234630` (`../launch.sh`, PUBLISH=1, 1 h 24 min);
  flank scan stages 1-7 and the report rebuild afterwards (`../fs_report.sh`), as for `rsi_main_20261016`.

## Verdicts (every copy; `loci_verdict.tsv`, per subfamily `flank_verdict.tsv`)
Owner's rule (2026-10-18): **SINE** = flanks tested unique (flankscan stage 9); **grey zone** = flanks masked as repeats or
untestable, or the whole subfamily could not be tested; **not SINE** = flanks shared with another copy (twin, array), the
same flank at the same place in >= 10 copies, or the subfamily's flanks still shared at the expansion limit.

| | copies |
|---|---|
| assigned | 63,769 |
| SINE loci | 51,989 |
| grey zone | 9,199 |
| not SINE (assigned) | 2,581 |

Whole subfamilies in the grey zone: MEG-RS (730; 3' side: no flank left to scan at the expansion limit), and MEG-RL (9),
MEG-TR (16), r2_3seqs (18), r4_32seqs (12): too few copies for the expansion test.

The top100 / rand100 / subfamily plates in `alignments/` are drawn from the SINE loci only: all 2,872 copy rows of the
top100 and rand100 plates have verdict SINE (checked against `loci_verdict.tsv`). Subfamilies without SINE loci have no
plates. Copy numbers in the report tables still count assigned copies; the report's overview and the per-subfamily table
give the SINE / grey / not-SINE split.

## Contents
| Path | What |
|---|---|
| `report.html` | full report (overview with the flank verdict, assignment, divergence, plates, element hierarchy, length versions, consensus audit) |
| `loci_verdict.tsv`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `alignments/` | top 100, random 100 and subfamily plates (SINE loci only) |
| `flankscan/` | hierarchy, peaks, candidates, chains, accepted candidates |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers (assigned copies) |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | stage-9 flank twins per family, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its own copies (computed from all assigned copies, not yet from SINE loci only) |
| `consensuses.clean.fa`, `manifest.txt` | the bank; tool versions and the command |

Not yet in this run: the unit-vs-composite assignment change (branch `coverage-tiebreak`, unmerged).

## What happened to r2 and r4 (checked 2026-10-18 against the peel's own loci)
The peel (2026-09-26) was made on chunks of a run that searched with **Rhin-1 only** (therioserver
`~/rhin/rsi/run_20260926_142634`; chunk loci in `genome.clean_step1/subfam_input/input_NNN.bnk`). Its loci are ~150 bp
monomer-length pieces. r2 = chunks 008, 312, 313 (150 loci); r4 = 32 chunks (1,600 loci). In this run 147 / 1,595 of them
are found (bedtools, >= 50 % overlap), but as parts of **whole composite copies**: median copy length 618 bp (C11, 608 bp
consensus) and 482 bp (P1, 448 bp), not ~150 bp.
- r2 loci -> C11 105, P1 32, P2 3, P34 1 (plus 6 unassigned). r2's consensus = C11 positions 193-338 at 99.3 %, i.e. C11's
  second module.
- r4 loci -> P1 1,111, P34 195, C11 156, P2 46, r5 8, r3 4 (plus 75 unassigned). r4's consensus = P1 200-335 at 97.8 %,
  C11 336-471 at 97.1 %, P34 140-275 at 96.3 %, r3 1-135 at 97.0 %: the r3-like module inside the composites.
So r2 and r4 were module-level groups sampled from inside composite elements; once the composites are in the bank, those
loci are assigned to the composites. Copies of r2 / r4 length on their own remain: 18 and 12 (too few for the flank
expansion test, hence grey). Verdicts of the peel loci in this run: r2 120 SINE / 17 grey / 10 not SINE; r4 1,242 / 195 / 158.
