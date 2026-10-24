# mev_main_20261024

*Myotis evotis* (Vespertilionidae; Vespertilionoidea, Yangochiroptera), genome `GCA_056320405.1`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/mev_main_20261024/` (run 2026-10-24 12:36:54 to 2026-10-24 16:07:34; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 225 | 196 | 26 | 3 | 244 |
| VES | VES | 316,847 | 221,409 | 42,340 | 53,098 | 315804 |
| MEG-RL | MEG | 4 | 0 | 4 | 0 | 6 |
| MEG-RS | MEG | 87 | 40 | 26 | 21 | 109 |
| MEG-T2 | MEG | 123 | 105 | 14 | 4 | 126 |
| MEG-TR | MEG | 17 | 17 | 0 | 0 | 19 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 225 | 196 | 26 | 3 |
| VES | 316,847 | 221,409 | 42,340 | 53,098 |
| MEG | 231 | 162 | 44 | 25 |

Copies in the list 457,942; assigned 317,303; not SINE in total (assigned + not assigned) 193,765.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
