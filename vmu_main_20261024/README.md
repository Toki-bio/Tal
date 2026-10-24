# vmu_main_20261024

*Vespertilio murinus* (Vespertilionidae; Vespertilionoidea, Yangochiroptera), genome `GCA_963924515.2`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/vmu_main_20261024/` (run 2026-10-24 13:38:14 to 2026-10-24 18:11:23; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 96 | 82 | 14 | 0 | 107 |
| VES | VES | 421,304 | 292,778 | 65,428 | 63,098 | 431357 |
| MEG-RL | MEG | 8 | 8 | 0 | 0 | 6 |
| MEG-RS | MEG | 158 | 52 | 15 | 91 | 622 |
| MEG-T2 | MEG | 263 | 85 | 81 | 97 | 327 |
| MEG-TR | MEG | 22 | 1 | 0 | 21 | 19 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 96 | 82 | 14 | 0 |
| VES | 421,304 | 292,778 | 65,428 | 63,098 |
| MEG | 451 | 146 | 96 | 209 |

Copies in the list 563,515; assigned 421,851; not SINE in total (assigned + not assigned) 204,971.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
