# rre_main_20261022

*Rhinolophus rex* (rre, GCA_041825415.1): complete SINEderella run with the *R. sinicus* 21-consensus bank (`bank21.fa`, the bank of
`rsi_main_20261021b`) and the same family table (`rsi_families.tsv`: Rhin = the r# subfamilies and the composites built from
them, MEG = the four MEG-* consensuses). The previous run of this genome with that bank is `rre_cross` (2026-10-06, code from
before the flank verdicts).

- Code: SINEderella `main` b0fc01e (verdicts in every run, SINE-discriminator modules in `disc/`, packed plate flanks,
  divergence against the searched consensus, small families tested, page table by family).
- Server: therioserver `/home/toki/rhin/rre_main_20261022/run_20261022_111022` (`../launch.sh`, `--families`, PUBLISH=1,
  THREADS=32), then flank-scan stages 1-7 and the report rebuild through `publish_run.sh`.
- Assembly QC: SOLID. Plates: 0 copy rows with a gap inside a flank.

## Verdicts (every copy; `loci_verdict.tsv`, per subfamily `flank_verdict.tsv`; rules: SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 97,996 |
| assigned | 92,854 |
| SINE loci | 76,619 |
| grey zone | 10,897 |
| not SINE (assigned + not assigned) | 10,480 |

Per family (assigned copies):

| family | subfamilies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| Rhin | 17 | 92,704 | 76,479 | 10,892 | 5,333 |
| MEG | 4 | 150 | 140 | 5 | 5 |

Every subfamily's flanks passed the subfamily test (all 21 interspersed).

## Contents
Same layout as `../rsi_main_20261021b/`: `report.html`, `loci_verdict.tsv`, `flank_verdict.tsv`, `alignments/` (plates of SINE
loci, topped up with grey copies marked [grey]), `flankscan/`, `satellites/`, `length_variants/`, `consensus_audit/`,
`consensuses.clean.fa`, `manifest.txt`, `rsi_families.tsv`, `families.summary.tsv`.
