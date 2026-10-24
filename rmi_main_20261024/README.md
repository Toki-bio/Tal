# rmi_main_20261024

*Rhinopoma microphyllum* (Rhinopomatidae; Rhinolophoidea, Yinpterochiroptera), genome `GCA_043880545.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/rmi_main_20261024/` (run 2026-10-24 08:44:25 to 2026-10-24 10:44:09; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 289 | 29 | 255 | 5 | 291 |
| VES | VES | 27,717 | 11,122 | 2,674 | 13,921 | 27163 |
| MEG-RL | MEG | 17 | 15 | 2 | 0 | 4 |
| MEG-RS | MEG | 24 | 19 | 5 | 0 | 86 |
| MEG-T2 | MEG | 257 | 174 | 76 | 7 | 211 |
| MEG-TR | MEG | 26 | 20 | 4 | 2 | 23 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 289 | 29 | 255 | 5 |
| VES | 27,717 | 11,122 | 2,674 | 13,921 |
| MEG | 324 | 228 | 87 | 9 |

Copies in the list 28,378; assigned 28,330; not SINE in total (assigned + not assigned) 13,983.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
