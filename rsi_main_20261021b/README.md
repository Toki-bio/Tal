# rsi_main_20261021b

*Rhinolophus sinicus* (rsi), complete SINEderella rerun: same genome (`GCF_057876585.1`), same 21-consensus bank and
settings as `rsi_main_20261021`, now with a family table (`rsi_families.tsv`: the r# subfamilies and the composites built
from them = **Rhin**, the four MEG-* consensuses = **MEG**).

- Code: SINEderella `main` def2b8a (2026-10-21). New against f442603:
  - plates keep their flanks packed against the element. Until now the last publish step (`tools/add_seed_row.py`)
    relabelled aligned element columns as flank where the copies reach past the searched consensus, so flanks showed gap
    columns (r10_19seqs, r6, r7, MEG-RS: every copy). This run: 0 of 3,258 plate copy rows with a gap inside a flank;
  - divergence plots measure each copy against the consensus it was searched with. Until now a bank rebuilt from the
    copies was used, which could differ: r1_r2_r3_r3_C11 showed 23 % divergence (median identity 76.5 %), against the
    searched consensus it is 6.5 % (93.5 %); r1_r3_P1 80.1 % -> 94.1 %;
  - wider page; the subfamily table has consensus length, assigned, SINE loci, grey and not-SINE columns, grouped
    under one row per family with the family's totals.
- Server: therioserver `/home/toki/rhin/rsi_main_20261021b/run_20261021_192808` (`../launch.sh`, `--families`, PUBLISH=1,
  THREADS=32; run 19:28-21:56, flank scan stages 1-7 and the report rebuild until 22:16).
- Assembly QC: SOLID.

## Verdicts (every copy; `loci_verdict.tsv`, per subfamily `flank_verdict.tsv`)
**SINE** = flanks tested unique; **grey zone** = flanks masked or untestable, or the subfamily could not be tested;
**not SINE** = flanks shared with another copy, the same flank at the same place in >= 10 copies, the subfamily's flanks
still shared at the expansion limit, or not assigned. Rules: SINEderella `docs/VERDICTS.md`.

| | rsi_main_20261021 | this run |
|---|---|---|
| copies in the list | 67,110 | 67,113 |
| assigned | 63,664 | 63,816 |
| SINE loci | 51,978 | 52,100 |
| grey zone | 9,157 | 9,173 |
| not SINE (assigned + not assigned) | 5,975 | 5,840 |

Per family (assigned copies: SINE / grey / not SINE):

| family | subfamilies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| Rhin | 17 | 62,950 | 51,980 | 8,438 | 2,532 |
| MEG | 4 | 866 | 120 | 735 | 11 |

The only subfamily grey as a whole is MEG-RS (733): its 3' side has no flank left to scan at the expansion limit (the
tandem-array family). The differences from `rsi_main_20261021` are within the run-to-run spread of the assignment vote
(ssearch36 `-z 11`).

## Contents
| Path | What |
|---|---|
| `report.html` | full report (subfamily table by family, plates, element hierarchy, consensus similarity, tRNA / 5S / 7SL origin, divergence, step1-2 tables) |
| `loci_verdict.tsv`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `rsi_families.tsv`, `families.summary.tsv` | the family table given to the run, and the family roll-up summary |
| `alignments/` | top 100, random 100 and subfamily plates (SINE loci, topped up with grey copies marked [grey]) |
| `flankscan/` | hierarchy, peaks, candidates, chains, accepted candidates |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers (assigned copies) |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | stage-9 flank twins per family, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its own copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composite consensuses read as units; the bank drawn as elementary blocks |
| `consensuses.clean.fa`, `manifest.txt` | the bank; tool versions and the command |

The r9 / r7 / r5 / r8 / r6 + r1 / r5 heads + tRNA-Ile-AAT alignment is consensus-level and stays in
`../rsi_main_20261018/alignments/rsi_r_units_heads_tRNA-Ile.aln.fa`.
