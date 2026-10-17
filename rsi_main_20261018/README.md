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
