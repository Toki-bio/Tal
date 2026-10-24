# tpe_main_20261024

*Triaenops persicus* (Rhinonycteridae; Rhinolophoidea, Yinpterochiroptera), genome `GCA_985607305.1`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/tpe_main_20261024/` (run 2026-10-24 12:04:42 to 2026-10-24 12:58:13; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 31,312 | 23,179 | 3,766 | 4,367 | 31409 |
| VES | VES | 31 | 25 | 6 | 0 | 36 |
| MEG-RL | MEG | 7 | 6 | 1 | 0 | 16 |
| MEG-RS | MEG | 25 | 18 | 0 | 7 | 23 |
| MEG-T2 | MEG | 93 | 91 | 2 | 0 | 104 |
| MEG-TR | MEG | 22 | 19 | 0 | 3 | 20 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 31,312 | 23,179 | 3,766 | 4,367 |
| VES | 31 | 25 | 6 | 0 |
| MEG | 147 | 134 | 3 | 10 |

Copies in the list 55,466; assigned 31,490; not SINE in total (assigned + not assigned) 28,353.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
