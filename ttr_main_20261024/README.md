# ttr_main_20261024

*Thyroptera tricolor* (Thyropteridae; Noctilionoidea, Yangochiroptera), genome `GCA_985616965.1`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/ttr_main_20261024/` (run 2026-10-24 11:04:42 to 2026-10-25 04:09:45; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 128 | 97 | 29 | 2 | 177 |
| VES | VES | 866,561 | 567,602 | 149,709 | 149,250 | 868458 |
| MEG-RL | MEG | 5 | 0 | 5 | 0 | 1 |
| MEG-RS | MEG | 38 | 19 | 13 | 6 | 45 |
| MEG-T2 | MEG | 65 | 60 | 3 | 2 | 88 |
| MEG-TR | MEG | 16 | 14 | 2 | 0 | 26 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 128 | 97 | 29 | 2 |
| VES | 866,561 | 567,602 | 149,709 | 149,250 |
| MEG | 124 | 93 | 23 | 8 |

Copies in the list 961,224; assigned 866,813; not SINE in total (assigned + not assigned) 243,671.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
