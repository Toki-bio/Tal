# cse_main_20261024

*Cistugo seabrae* (Cistugidae; Vespertilionoidea, Yangochiroptera), genome `GCA_061228675.1`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/cse_main_20261024/` (run 2026-10-24 15:54:45 to 2026-10-24 20:55:11; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 239 | 197 | 42 | 0 | 270 |
| VES | VES | 570,370 | 421,455 | 77,293 | 71,622 | 568685 |
| MEG-RL | MEG | 8 | 8 | 0 | 0 | 6 |
| MEG-RS | MEG | 36 | 32 | 4 | 0 | 119 |
| MEG-T2 | MEG | 115 | 106 | 5 | 4 | 207 |
| MEG-TR | MEG | 12 | 11 | 1 | 0 | 14 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 239 | 197 | 42 | 0 |
| VES | 570,370 | 421,455 | 77,293 | 71,622 |
| MEG | 171 | 157 | 10 | 4 |

Copies in the list 626,239; assigned 570,780; not SINE in total (assigned + not assigned) 127,085.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
