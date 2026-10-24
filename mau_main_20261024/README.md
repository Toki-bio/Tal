# mau_main_20261024

*Myzopoda aurita* (Myzopodidae; Emballonuroidea, Yangochiroptera), genome `GCA_985615965.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/mau_main_20261024/` (run 2026-10-24 08:44:25 to 2026-10-24 12:36:53; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 176 | 160 | 14 | 2 | 196 |
| VES | VES | 425,267 | 325,755 | 53,005 | 46,507 | 424776 |
| MEG-RL | MEG | 13 | 12 | 1 | 0 | 13 |
| MEG-RS | MEG | 76 | 65 | 11 | 0 | 97 |
| MEG-T2 | MEG | 25 | 24 | 1 | 0 | 40 |
| MEG-TR | MEG | 12 | 10 | 0 | 2 | 17 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 176 | 160 | 14 | 2 |
| VES | 425,267 | 325,755 | 53,005 | 46,507 |
| MEG | 126 | 111 | 13 | 2 |

Copies in the list 446,824; assigned 425,569; not SINE in total (assigned + not assigned) 67,766.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
