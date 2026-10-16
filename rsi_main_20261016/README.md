# rsi_main_20261016

*Rhinolophus sinicus* (rsi), complete SINEderella run on the current pipeline, same 21-consensus bank as `rsi_final`
(r1-r10, MEG-RL/RS/T2/TR, seven composites; `consensuses.clean.fa`).

- Code: SINEderella `main` 296fa33 (2026-10-14 audit fixes, strand-aware merge, kind-B array check with short pieces,
  satellite screen on, SubFam 1.5.0). Assembly `GCF_057876585.1`. THREADS 32, run 1 h 26 min.
- Server: therioserver `/home/toki/rhin/rsi_main_20261016/run_20261016_194647` (launcher `../launch.sh`).
- Not in this run (they were added to `rsi_final` by hand): the separate flank scan (stages 1-7) and the stage-8 plates.
- The unit-vs-composite fix (branch `coverage-tiebreak`) is NOT in this run; it is not ready (see
  `rsi_composites_test/`).

## Contents
| Path | What |
|---|---|
| `report.html` | full report (assignment, divergence, plates linking to `alignments/`, length versions, consensus audit) |
| `alignments/` | top 100, random 100 and subfamily plates (60 files) |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |
| `satellites/indication.tsv`, `satellites/loci.bed` | satellite screen: SINE-derived tandem loci (kind A) and verified arrays (kind B) per consensus |
| `array_flag.tsv`, `flank_twins.tsv` | tandem-array share per family; copies with shared flanks (flankscan fs9) |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its own copies |
| `consensuses.clean.fa`, `manifest.txt` | the bank; tool versions and the command |

## Against rsi_final (assigned copies)
Every family within ~1 % except:
- **MEG-RS 1,715 -> 730**: the satellite screen (not run in `rsi_final`) removed 15 verified MEG-RS tandem arrays
  (`satellites/indication.tsv`). This is the intended effect: `rsi_final` showed 1,592 of 1,715 MEG-RS copies in arrays.
- **r6 6,585 -> 7,019 (+434)**: r6 is the family whose copies flip most between runs of the voting step (294 copies
  firm -> not firm between two step-2 runs on identical input, 2026-10-16); the strand-aware merge also changes loci.
- Total assigned 64,555 -> 63,997.
