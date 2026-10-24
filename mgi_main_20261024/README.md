# mgi_main_20261024

*Macroderma gigas* (Megadermatidae; Rhinolophoidea, Yinpterochiroptera), genome `GCA_055471755.1`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/mgi_main_20261024/` (run 2026-10-24 13:24:40 to 2026-10-24 13:38:08; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 292 | 268 | 22 | 2 | 303 |
| VES | VES | 12 | 11 | 1 | 0 | 12 |
| MEG-RL | MEG | 5 | 0 | 5 | 0 | 7 |
| MEG-RS | MEG | 9 | 9 | 0 | 0 | 9 |
| MEG-T2 | MEG | 39 | 39 | 0 | 0 | 70 |
| MEG-TR | MEG | 12 | 12 | 0 | 0 | 16 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 292 | 268 | 22 | 2 |
| VES | 12 | 11 | 1 | 0 |
| MEG | 65 | 60 | 5 | 0 |

Copies in the list 486; assigned 369; not SINE in total (assigned + not assigned) 119.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
