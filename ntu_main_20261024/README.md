# ntu_main_20261024

*Natalus tumidirostris* (Natalidae; Vespertilionoidea, Yangochiroptera), genome `GCA_985614685.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/ntu_main_20261024/` (run 2026-10-24 10:24:23 to 2026-10-24 13:57:10; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 649 | 578 | 71 | 0 | 734 |
| VES | VES | 344,170 | 221,847 | 71,847 | 50,476 | 345052 |
| MEG-RL | MEG | 9 | 9 | 0 | 0 | 6 |
| MEG-RS | MEG | 37 | 26 | 2 | 9 | 50 |
| MEG-T2 | MEG | 71 | 67 | 0 | 4 | 77 |
| MEG-TR | MEG | 11 | 9 | 2 | 0 | 14 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 649 | 578 | 71 | 0 |
| VES | 344,170 | 221,847 | 71,847 | 50,476 |
| MEG | 128 | 111 | 4 | 13 |

Copies in the list 367,068; assigned 344,947; not SINE in total (assigned + not assigned) 72,610.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
