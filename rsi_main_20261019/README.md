# rsi_main_20261019

*Rhinolophus sinicus* (rsi), complete SINEderella rerun: same genome (`GCF_057876585.1`), same 21-consensus bank and
settings as `rsi_main_20261018`; only the code changed.

- Code: SINEderella `main` 7cb38dc (2026-10-19). New against 6411da2: the GLM-audit fixes (docs/AUDIT_FIXES_6.md and
  batches 3-5: failed ssearch36 calls stop the run instead of reading as "no hits", satellite stage fixes, plate verdict
  look-up by the row's own subfamily, edge-profile fixes), the strict consensus rebuild (a rebuild that did not converge
  counts as failed), `--add` searches on the source run's mask, and the assembly-QC line on the report.
- Server: therioserver `/home/toki/rhin/rsi_main_20261019/run_20261019_151118` (`../launch.sh`, PUBLISH=1, THREADS=32;
  run 15:11-18:33, flank scan stages 1-7 and the report rebuild through `publish_run.sh` until 18:57).
- Assembly QC: SOLID.

## Verdicts (every copy; `loci_verdict.tsv`, per subfamily `flank_verdict.tsv`)
Owner's rule (2026-10-18): **SINE** = flanks tested unique (flankscan stage 9); **grey zone** = flanks masked as repeats or
untestable, or the whole subfamily could not be tested; **not SINE** = flanks shared with another copy (twin, array), the
same flank at the same place in >= 10 copies, or the subfamily's flanks still shared at the expansion limit.

| | rsi_main_20261018 | this run |
|---|---|---|
| copies in the list | 67,118 | 67,120 |
| assigned | 63,769 | 63,837 |
| SINE loci | 51,989 | 52,075 |
| grey zone | 9,199 | 9,228 |
| not SINE (assigned) | 2,581 | 2,534 |
| not assigned (no unanimous vote / low bitscore) | 3,349 | 3,283 |

Per subfamily the changes are within the run-to-run spread of the assignment vote (ssearch36 `-z 11`), except MEG-TR:
20 assigned copies here against 16 before. Flankscan stage 9 tests families with at least 20 assigned copies
(`flankscan/fs9_run.sh`, MIN_COPIES=20), so MEG-TR was tested this time: 18 SINE, 2 not SINE (flank twins); in
`rsi_main_20261018` all 16 were grey. The remaining whole-subfamily grey zones: MEG-RS (730), and MEG-RL (13),
r2_3seqs (18), r4_32seqs (13), which have fewer than 20 assigned copies.

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
