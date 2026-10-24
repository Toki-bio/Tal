# tbr_main_20261024

*Tadarida brasiliensis* (Molossidae; Vespertilionoidea, Yangochiroptera), genome `GCA_030848825.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/tbr_main_20261024/` (run 2026-10-24 09:27:15 to 2026-10-24 17:29:07; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 366 | 306 | 58 | 2 | 437 |
| VES | VES | 620,807 | 373,297 | 136,235 | 111,275 | 621128 |
| MEG-RL | MEG | 11 | 10 | 1 | 0 | 9 |
| MEG-RS | MEG | 53 | 37 | 5 | 11 | 249 |
| MEG-T2 | MEG | 78 | 76 | 2 | 0 | 92 |
| MEG-TR | MEG | 21 | 14 | 1 | 6 | 24 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 366 | 306 | 58 | 2 |
| VES | 620,807 | 373,297 | 136,235 | 111,275 |
| MEG | 163 | 137 | 9 | 17 |

Copies in the list 639,492; assigned 621,336; not SINE in total (assigned + not assigned) 129,450.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
