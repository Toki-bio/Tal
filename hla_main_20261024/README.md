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
