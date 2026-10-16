# rsi_main_20261016

*Rhinolophus sinicus* (rsi), complete SINEderella run on the current pipeline, same 21-consensus bank as `rsi_final`
(r1-r10, MEG-RL/RS/T2/TR, seven composites; `consensuses.clean.fa`).

- Code: SINEderella `main` 296fa33 (2026-10-14 audit fixes, strand-aware merge, kind-B array check with short pieces,
  satellite screen on, SubFam 1.5.0). Assembly `GCF_057876585.1`. THREADS 32, run 1 h 26 min.
- Server: therioserver `/home/toki/rhin/rsi_main_20261016/run_20261016_194647` (launcher `../launch.sh`).
- Flank scan stages 1-7 (`flankscan/fs_all.sh`, 18 min) were run on this run afterwards, as for `rsi_final`, and the report was
  rebuilt with them: it has the **Element hierarchy** section again (`flankscan/hierarchy.tsv`, 14 elements; `rsi_final` 13).
  Stage-8 plates were not made for this run.
- The unit-vs-composite fix (branch `coverage-tiebreak`) is NOT in this run; it is not ready (see
  `rsi_composites_test/`).

## Contents
| Path | What |
|---|---|
| `composites.html` | scheme of the 7 composites (parts at their placed positions) and the units, with tracks for stretches similar to differently named units, and the alignment of all 21 consensuses (draft design, 2026-10-17) |
| `consensuses.aln.fa` | MAFFT L-INS-i alignment of the 21 consensuses |
| `report.html` | full report (assignment, divergence, plates linking to `alignments/`, element hierarchy, length versions, consensus audit) |
| `flankscan/` | hierarchy, peaks, candidates, re-assignment and chain tables, accepted candidates (`cand/` alignments) |
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
