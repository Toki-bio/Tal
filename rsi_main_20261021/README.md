# rsi_main_20261021

*Rhinolophus sinicus* (rsi), complete SINEderella rerun: same genome (`GCF_057876585.1`), same 21-consensus bank and
settings as `rsi_main_20261019`; only the code changed.

- Code: SINEderella `main` f442603 (2026-10-21). New against 7cb38dc: small families are flank-tested. The subfamily
  test (step7) starts at 7 copies instead of 20 and compares every pair of copies in small families; the per-copy flank-twin
  test (flankscan stage 9) starts at 2 copies instead of 20. Until now families under 20 copies were grey only for their
  size, and a family near 20 could change verdict between runs (MEG-TR: 16 copies grey in `rsi_main_20261018`, 20 tested
  in `rsi_main_20261019`). Test: SINEderella `tests/test_step7_small_families.sh`.
- Server: therioserver `/home/toki/rhin/rsi_main_20261021/run_20261021_135028` (`../launch.sh`, PUBLISH=1, THREADS=32;
  run 13:50-15:29, flank scan stages 1-7 and the report rebuild through `publish_run.sh` until 15:50).
- Assembly QC: SOLID.

## Verdicts (every copy; `loci_verdict.tsv`, per subfamily `flank_verdict.tsv`)
Owner's rule (2026-10-18): **SINE** = flanks tested unique (flankscan stage 9); **grey zone** = flanks masked as repeats or
untestable, or the whole subfamily could not be tested; **not SINE** = flanks shared with another copy (twin, array), the
same flank at the same place in >= 10 copies, or the subfamily's flanks still shared at the expansion limit.

| | rsi_main_20261019 | this run |
|---|---|---|
| copies in the list | 67,120 | 67,110 |
| assigned | 63,837 | 63,664 |
| SINE loci | 52,075 | 51,978 |
| grey zone | 9,228 | 9,157 |
| not SINE (assigned) | 2,534 | 2,529 |
| not assigned (no unanimous vote / low bitscore) | 3,283 | 3,446 |

The small families, tested for the first time (SINE / grey / not SINE of the assigned copies):

| subfamily | assigned | rsi_main_20261019 | this run |
|---|---|---|---|
| MEG-RL | 14 | 0 / 13 / 0 (untested) | 10 / 1 / 3 |
| MEG-TR | 18 | 18 / 0 / 2 | 16 / 0 / 2 |
| r2_3seqs | 18 | 0 / 18 / 0 (untested) | 17 / 1 / 0 |
| r4_32seqs | 14 | 0 / 13 / 0 (untested) | 12 / 2 / 0 |

(The 2026-10-19 column counts that run's assigned copies: MEG-RL 13, MEG-TR 20, r2 18, r4 13.) The only subfamily still grey
as a whole is MEG-RS (732): its 3' side has no flank left to scan at the expansion limit (the tandem-array family). The other
differences between the two runs are within the run-to-run spread of the assignment vote (ssearch36 `-z 11`).

## Contents
| Path | What |
|---|---|
| `report.html` | full report (overview with the flank verdict and the assembly-QC line, plates, element hierarchy, consensus similarity, tRNA / 5S / 7SL origin, divergence, step1-2 tables) |
| `loci_verdict.tsv`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `alignments/` | top 100, random 100 and subfamily plates (SINE loci only) |
| `flankscan/` | hierarchy, peaks, candidates, chains, accepted candidates |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers (assigned copies) |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | stage-9 flank twins per family, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its own copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composite consensuses read as units; the bank drawn as elementary blocks |
| `consensuses.clean.fa`, `manifest.txt` | the bank; tool versions and the command |

The r9 / r7 / r5 / r8 / r6 + r1 / r5 heads + tRNA-Ile-AAT alignment is consensus-level and stays in
`../rsi_main_20261018/alignments/rsi_r_units_heads_tRNA-Ile.aln.fa`.
