# lly_main_20261024

*Lyroderma lyra* (Megadermatidae; Rhinolophoidea, Yinpterochiroptera), genome `GCA_004026885.1`. Complete SINEderella rerun on
main 8d4961c: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/lly_main_20261024/` (run 2026-10-24 17:48:50 to 2026-10-24 18:16:48; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 319 | 287 | 28 | 4 | 290 |
| VES | VES | 10 | 9 | 1 | 0 | 11 |
| MEG-RL | MEG | 243 | 42 | 84 | 117 | 23 |
| MEG-RS | MEG | 21 | 16 | 5 | 0 | 17 |
| MEG-T2 | MEG | 82 | 75 | 5 | 2 | 85 |
| MEG-TR | MEG | 24 | 23 | 1 | 0 | 24 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 319 | 287 | 28 | 4 |
| VES | 10 | 9 | 1 | 0 |
| MEG | 370 | 156 | 95 | 119 |

Copies in the list 1,078; assigned 699; not SINE in total (assigned + not assigned) 502.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
