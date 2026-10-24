# hla_main_20261024

*Hipposideros larvatus* (Hipposideridae; Rhinolophoidea, Yinpterochiroptera), genome `GCA_031876335.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/hla_main_20261024/` (run 2026-10-24 08:44:25 to 2026-10-24 09:27:14; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 20,401 | 13,472 | 3,245 | 3,684 | 20602 |
| VES | VES | 20 | 17 | 3 | 0 | 22 |
| MEG-RL | MEG | 423 | 0 | 0 | 423 | 4029 |
| MEG-RS | MEG | 30 | 19 | 3 | 8 | 22 |
| MEG-T2 | MEG | 80 | 79 | 1 | 0 | 90 |
| MEG-TR | MEG | 16 | 10 | 2 | 4 | 527 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 20,401 | 13,472 | 3,245 | 3,684 |
| VES | 20 | 17 | 3 | 0 |
| MEG | 549 | 108 | 6 | 435 |

Copies in the list 71,135; assigned 20,970; not SINE in total (assigned + not assigned) 54,284.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).

## Why MEG-RL and MEG-TR dropped against the 2026-09-27 run
The satellite screen (added to SINEderella after that run) removed 5,955 hits lying in 838 verified tandem arrays
(`satellites/indication.tsv`: 89.5 % of MEG-RS hits and 88.9 % of MEG-TR hits in arrays; arrays 233 bp to 116 kb, median 38 kb).
Of the 4,029 firmly assigned MEG-RL copies of the old run, 3,572 lie inside these arrays; of 527 MEG-TR copies, 511.
Example, CM063088.1:109.71-109.80 Mb: 14 old MEG-RL copies at a regular spacing of 5.42-5.46 kb (gaps of about 10.9 and
13.3 kb where a unit is missing), i.e. one MEG copy per 5.4 kb array unit. The 423 MEG-RL copies left are verdicted not a
SINE at subfamily level: their flanks are still shared at the expansion limit. Rhin-1 is unchanged (20,602 -> 20,401).
