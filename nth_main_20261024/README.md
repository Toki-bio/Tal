# nth_main_20261024

*Nycteris thebaica* (Nycteridae; Emballonuroidea, Yangochiroptera), genome `GCA_985607365.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/nth_main_20261024/` (run 2026-10-24 08:44:25 to 2026-10-24 09:30:54; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 1,180 | 865 | 286 | 29 | 1283 |
| VES | VES | 23,742 | 15,655 | 2,725 | 5,362 | 24084 |
| MEG-RL | MEG | 4 | 0 | 4 | 0 | 4 |
| MEG-RS | MEG | 62 | 33 | 21 | 8 | 7 |
| MEG-T2 | MEG | 103 | 92 | 11 | 0 | 98 |
| MEG-TR | MEG | 19 | 16 | 3 | 0 | 23 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 1,180 | 865 | 286 | 29 |
| VES | 23,742 | 15,655 | 2,725 | 5,362 |
| MEG | 188 | 141 | 39 | 8 |

Copies in the list 25,300; assigned 25,110; not SINE in total (assigned + not assigned) 5,589.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
